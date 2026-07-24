// Named parameter presets + last-used state, persisted in userData (GUI P2).
//
// One small JSON file. Writes are atomic (write a temp, rename over): a crash
// mid-write leaves the previous good file, never a half-written one. A corrupt
// or absent file yields empty defaults rather than throwing -- losing presets is
// bad, but crashing the app on launch because of them is worse.

import { app } from 'electron'
import { existsSync, readFileSync, writeFileSync, renameSync, mkdirSync } from 'node:fs'
import { join, dirname } from 'node:path'

export interface Settings {
  schemaVersion: number
  /// The form state to restore on next launch (values keyed by CLI option name).
  lastUsed: Record<string, unknown> | null
  /// name -> parameter values.
  presets: Record<string, Record<string, unknown>>
}

const FILE = () => join(app.getPath('userData'), 'fastag-settings.json')
const EMPTY: Settings = { schemaVersion: 1, lastUsed: null, presets: {} }

export function loadSettings(): Settings {
  try {
    const p = FILE()
    if (!existsSync(p)) return { ...EMPTY }
    const raw = JSON.parse(readFileSync(p, 'utf8'))
    // Tolerate an older/partial shape rather than trust the file blindly.
    return {
      schemaVersion: 1,
      lastUsed: raw && typeof raw.lastUsed === 'object' ? raw.lastUsed : null,
      presets: raw && typeof raw.presets === 'object' && raw.presets ? raw.presets : {}
    }
  } catch {
    return { ...EMPTY }
  }
}

function save(s: Settings): boolean {
  try {
    const p = FILE()
    mkdirSync(dirname(p), { recursive: true })
    const tmp = `${p}.tmp`
    writeFileSync(tmp, JSON.stringify(s, null, 2), 'utf8')
    renameSync(tmp, p)
    return true
  } catch {
    return false
  }
}

export function saveLastUsed(values: Record<string, unknown>): boolean {
  const s = loadSettings()
  s.lastUsed = values
  return save(s)
}

export function savePreset(name: string, values: Record<string, unknown>): boolean {
  const clean = name.trim()
  if (!clean) return false
  const s = loadSettings()
  s.presets[clean] = values
  return save(s)
}

export function deletePreset(name: string): boolean {
  const s = loadSettings()
  if (!(name in s.presets)) return false
  delete s.presets[name]
  return save(s)
}
