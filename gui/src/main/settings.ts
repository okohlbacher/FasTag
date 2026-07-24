// Named parameter presets + last-used state, persisted in userData (GUI P2).
//
// One small JSON file. Writes are atomic (write a temp, rename over): a crash
// mid-write leaves the previous good file, never a half-written one. A corrupt
// or absent file yields empty defaults rather than throwing -- losing presets is
// bad, but crashing the app on launch because of them is worse.

import { app } from 'electron'
import { existsSync, readFileSync, writeFileSync, renameSync, mkdirSync } from 'node:fs'
import { join, dirname } from 'node:path'

// __proto__/constructor/prototype as a preset name or last-used key would poison
// Object.prototype when merged in the renderer. Keep only own, safe keys.
function sanitize<T extends object>(o: T): T {
  const out = Object.create(null) as Record<string, unknown>
  for (const [k, v] of Object.entries(o)) {
    if (k === '__proto__' || k === 'constructor' || k === 'prototype') continue
    out[k] = v
  }
  return out as T
}

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
  const p = FILE()
  try {
    if (!existsSync(p)) return { ...EMPTY }
    const raw = JSON.parse(readFileSync(p, 'utf8'))
    // Tolerate an older/partial shape rather than trust the file blindly.
    const presetsIn =
      raw && typeof raw.presets === 'object' && raw.presets ? sanitize(raw.presets as Record<string, Record<string, unknown>>) : {}
    const presets: Record<string, Record<string, unknown>> = {}
    for (const [name, values] of Object.entries(presetsIn)) {
      if (values && typeof values === 'object') presets[name] = sanitize(values)
    }
    return {
      schemaVersion: 1,
      lastUsed: raw && typeof raw.lastUsed === 'object' && raw.lastUsed ? sanitize(raw.lastUsed) : null,
      presets
    }
  } catch {
    // The file exists but didn't parse. Move it aside ONCE so the next save()
    // doesn't silently overwrite (and destroy) recoverable presets; then start
    // empty rather than crashing the app on launch.
    try {
      if (existsSync(p) && !existsSync(`${p}.corrupt`)) renameSync(p, `${p}.corrupt`)
    } catch {
      /* best effort -- a failed backup must not itself throw */
    }
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
  s.lastUsed = sanitize(values)
  return save(s)
}

export function savePreset(name: string, values: Record<string, unknown>): boolean {
  const clean = name.trim()
  if (!clean || clean === '__proto__' || clean === 'constructor' || clean === 'prototype') return false
  const s = loadSettings()
  s.presets[clean] = sanitize(values)
  return save(s)
}

export function deletePreset(name: string): boolean {
  const s = loadSettings()
  if (!(name in s.presets)) return false
  delete s.presets[name]
  return save(s)
}
