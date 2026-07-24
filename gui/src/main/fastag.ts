// FasTag CLI integration: locate the bundled (or dev) binary, run it safely,
// stream its stderr log, and cancel a run by killing its process tree.
//
// Security: the binary is spawned with shell:false and an explicit argv array —
// no shell string is ever built from user input, so there is nothing to inject.
// Every option value is passed as its own argv element.

import { spawn, spawnSync, ChildProcess } from 'node:child_process'
import { existsSync } from 'node:fs'
import { join, dirname } from 'node:path'
import { app } from 'electron'

export interface BinaryInfo {
  bin: string
  dataPath?: string
  source: 'env' | 'bundled' | 'path'
  ok: boolean
  version?: string
  detail: string
}

// The FasTag options the P0/P1 skeleton drives. This is the curated subset;
// the full self-describing param manifest (from -write_ini) comes in a later
// phase. Keeping it explicit means we never interpolate an unknown flag.
export interface RunParams {
  in: string
  out: string
  tag_length?: number
  fragment_tolerance?: number
  fragment_tolerance_unit?: 'ppm' | 'Da'
  max_tags?: number
  threads?: number
  progress?: boolean
}

const exeName = process.platform === 'win32' ? 'FasTag.exe' : 'FasTag'

// Resolve the binary in priority order:
//   1. FASTAG_BIN env override (developer pointing at a local build)
//   2. bundled inside the app under resources/fastag/bin (production)
//   3. bare name on PATH (last resort)
// The bundled layout mirrors the release tarball: bin/FasTag + share/OpenMS.
export function resolveBinary(): { bin: string; dataPath?: string; source: BinaryInfo['source'] } {
  const envBin = process.env.FASTAG_BIN
  if (envBin && existsSync(envBin)) {
    return { bin: envBin, dataPath: process.env.OPENMS_DATA_PATH, source: 'env' }
  }

  const root = app.isPackaged
    ? join(process.resourcesPath, 'fastag')
    : join(app.getAppPath(), 'resources', 'fastag')
  const bundled = join(root, 'bin', exeName)
  if (existsSync(bundled)) {
    const dataPath = join(root, 'share', 'OpenMS')
    return { bin: bundled, dataPath: existsSync(dataPath) ? dataPath : undefined, source: 'bundled' }
  }

  return { bin: exeName, source: 'path' }
}

// The architecture canary (P0): prove the resolved binary actually runs and
// links its libraries on THIS machine, before the UI offers to run anything.
// TOPP tools print their version banner on --help; capture it as proof of life.
export function probeBinary(): BinaryInfo {
  const { bin, dataPath, source } = resolveBinary()
  const env = { ...process.env }
  if (dataPath) env.OPENMS_DATA_PATH = dataPath

  try {
    const r = spawnSync(bin, ['--help'], { env, encoding: 'utf8', timeout: 15000 })
    if (r.error) {
      return { bin, dataPath, source, ok: false, detail: `cannot execute: ${r.error.message}` }
    }
    const text = `${r.stdout || ''}\n${r.stderr || ''}`
    // TOPP --help banner looks like: "FasTag -- <desc>" then "Version: X.Y.Z ..."
    const m = text.match(/Version:\s*([0-9][^\s,]*)/i)
    return {
      bin,
      dataPath,
      source,
      ok: true,
      version: m ? m[1] : undefined,
      detail: m ? `FasTag ${m[1]} (${source})` : `runs, version unknown (${source})`
    }
  } catch (e) {
    return { bin, dataPath, source, ok: false, detail: `probe failed: ${(e as Error).message}` }
  }
}

// Build the argv array from typed params. Only known flags are emitted, and
// each value is its own array element — never concatenated into a string.
function buildArgs(p: RunParams): string[] {
  const args = ['-in', p.in, '-out', p.out]
  if (p.tag_length != null) args.push('-tag_length', String(p.tag_length))
  if (p.fragment_tolerance != null) args.push('-fragment_tolerance', String(p.fragment_tolerance))
  // The tool's tolerance is unit-sensitive (default ppm); a value without its
  // unit is meaningless. Always pass the unit alongside the value.
  if (p.fragment_tolerance_unit != null) args.push('-fragment_tolerance_unit', p.fragment_tolerance_unit)
  if (p.max_tags != null) args.push('-max_tags', String(p.max_tags))
  if (p.threads != null) args.push('-threads', String(p.threads))
  if (p.progress) args.push('-progress')
  return args
}

export interface RunHandle {
  child: ChildProcess
  args: string[]
}

export interface RunCallbacks {
  onLog: (line: string) => void
  onProgress: (done: number, total: number) => void
  onExit: (code: number | null, signal: NodeJS.Signals | null) => void
  onError: (message: string) => void
}

// The CLI emits `FASTAG_PROGRESS done=<n> total=<n>` on stderr when -progress is
// set. total is 0 for the streaming mzPeak path (count not known ahead of time).
const PROGRESS_RE = /^FASTAG_PROGRESS done=(\d+) total=(\d+)/

// Spawn a run. Returns the handle so the caller can cancel it. stderr is the
// CLI's progress channel; we split it into lines and forward each.
export function runFastag(params: RunParams, cb: RunCallbacks): RunHandle {
  const { bin, dataPath } = resolveBinary()
  const env = { ...process.env }
  if (dataPath) env.OPENMS_DATA_PATH = dataPath

  const args = buildArgs(params)
  const child = spawn(bin, args, { env, windowsHide: true })

  let buf = ''
  const routeLine = (line: string): void => {
    const m = PROGRESS_RE.exec(line)
    if (m) cb.onProgress(Number(m[1]), Number(m[2]))
    else cb.onLog(line)
  }
  const pump = (chunk: Buffer): void => {
    buf += chunk.toString('utf8')
    let nl: number
    while ((nl = buf.indexOf('\n')) >= 0) {
      routeLine(buf.slice(0, nl))
      buf = buf.slice(nl + 1)
    }
  }
  child.stderr?.on('data', pump)
  child.stdout?.on('data', pump)
  child.on('error', (e) => cb.onError(e.message))
  child.on('close', (code, signal) => {
    if (buf.length) cb.onLog(buf)
    cb.onExit(code, signal)
  })

  return { child, args }
}

// Cancel: kill the whole process tree (the CLI may spawn OpenMP workers).
export function cancelRun(handle: RunHandle): void {
  const pid = handle.child.pid
  if (pid == null) return
  if (process.platform === 'win32') {
    spawnSync('taskkill', ['/pid', String(pid), '/t', '/f'])
  } else {
    handle.child.kill('SIGTERM')
    // Escalate if it does not exit promptly.
    setTimeout(() => {
      if (!handle.child.killed) handle.child.kill('SIGKILL')
    }, 2000)
  }
}
