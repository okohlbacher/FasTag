// The single bridge between the React frontend and the Rust backend. It presents
// the exact `window.fastag` shape the Electron build exposed, so App.tsx and its
// run/batch orchestration port over unchanged: only the plumbing underneath is
// Tauri (invoke + events) instead of Electron IPC.

import { invoke } from '@tauri-apps/api/core'
import { listen, type UnlistenFn } from '@tauri-apps/api/event'
import { open, save } from '@tauri-apps/plugin-dialog'
import { openUrl } from '@tauri-apps/plugin-opener'
import type {
  BinaryInfo,
  FastagApi,
  Preview,
  RunResult,
  RunStarted,
  Settings,
  SpeciesReport,
  TaxdbInfo
} from './types'

const SPECTRA = [{ name: 'Spectra', extensions: ['mzML', 'mzml', 'mzpeak', 'mzPeak'] }]

// listen() is async but the renderer wants a synchronous unsubscribe (it returns
// one from a useEffect). Hand back a function that unlistens once the async
// subscription resolves, and no-ops safely if it fires before then.
function subscribe<T>(event: string, cb: (payload: T) => void): () => void {
  let un: UnlistenFn | null = null
  let cancelled = false
  listen<T>(event, (e) => cb(e.payload)).then((u) => {
    if (cancelled) u()
    else un = u
  })
  return () => {
    cancelled = true
    if (un) un()
  }
}

const api: FastagApi = {
  probe: () => invoke<BinaryInfo>('probe'),

  // Native dialogs live in the frontend via the dialog plugin; the backend only
  // ever receives already-chosen paths.
  pickInput: async () => {
    const r = await open({ multiple: false, filters: SPECTRA })
    return typeof r === 'string' ? r : null
  },
  pickInputs: async () => {
    const r = await open({ multiple: true, filters: SPECTRA })
    return Array.isArray(r) ? r : []
  },
  pickOutput: async (defaultPath?: string) => {
    const r = await save({ defaultPath, filters: [{ name: 'Tags (TSV)', extensions: ['tsv'] }] })
    return typeof r === 'string' ? r : null
  },

  run: (params: Record<string, unknown>) => invoke<RunStarted>('run', { params }),
  cancel: () => invoke<{ cancelled: boolean }>('cancel'),
  preview: (path: string, maxRows?: number) => invoke<Preview>('preview', { path, maxRows }),
  species: (path: string) => invoke<SpeciesReport | null>('species', { path }),
  taxdbInfo: (explicit?: string) => invoke<TaxdbInfo | null>('taxdb_info', { explicit }),

  loadSettings: () => invoke<Settings>('load_settings'),
  saveLast: (values: Record<string, unknown>) => invoke<boolean>('save_last', { values }),
  savePreset: (name: string, values: Record<string, unknown>) =>
    invoke<boolean>('save_preset', { name, values }),
  deletePreset: (name: string) => invoke<boolean>('delete_preset', { name }),

  // Built from a fixed template here (the renderer passes a number, never a URL),
  // so a crafted value cannot redirect the user elsewhere.
  openTaxon: async (taxid: number) => {
    const id = Number(taxid)
    if (!Number.isInteger(id) || id <= 0) return false
    await openUrl(`https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=${id}`)
    return true
  },

  onLog: (cb) => subscribe<string>('fastag:log', cb),
  onProgress: (cb) => subscribe<{ done: number; total: number }>('fastag:progress', cb),
  onDone: (cb) => subscribe<RunResult>('fastag:done', cb)
}

declare global {
  interface Window {
    fastag: FastagApi
  }
}

window.fastag = api
