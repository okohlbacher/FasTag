// Shared shapes between preload and renderer.

export interface BinaryInfo {
  bin: string
  dataPath?: string
  source: 'env' | 'bundled' | 'path'
  ok: boolean
  version?: string
  detail: string
}

export interface Preview {
  header: string[]
  rows: string[][]
  truncated: boolean
  shown: number
}

export interface RunResult {
  ok: boolean
  code: number | null
  message?: string
  // Set on results delivered via fastag:done; correlates the terminal event to
  // the run the renderer is awaiting. Absent on synthetic failure results the
  // renderer builds itself (e.g. a run that never started).
  runId?: number
}

export interface Taxon {
  rank: string
  taxid: number
  name: string
  observed: number
  expected: number
  enrichment: number
  logP: number
  q: number
}

export interface SpeciesReport {
  path: string
  taxa: Taxon[]
  empty: boolean
}

export interface TaxdbInfo {
  path: string
  k: number
  kmers: number
}

export interface Settings {
  schemaVersion: number
  lastUsed: Record<string, unknown> | null
  presets: Record<string, Record<string, unknown>>
}

export interface FastagApi {
  loadSettings: () => Promise<Settings>
  saveLast: (values: Record<string, unknown>) => Promise<boolean>
  savePreset: (name: string, values: Record<string, unknown>) => Promise<boolean>
  deletePreset: (name: string) => Promise<boolean>
  taxdbInfo: (explicit?: string) => Promise<TaxdbInfo | null>
  species: (path: string) => Promise<SpeciesReport | null>
  openTaxon: (taxid: number) => Promise<boolean>
  probe: () => Promise<BinaryInfo>
  pickInput: () => Promise<string | null>
  pickInputs: () => Promise<string[]>
  pickOutput: (defaultPath?: string) => Promise<string | null>
  run: (params: Record<string, unknown>) => Promise<{ started: boolean; reason?: string; runId?: number }>
  cancel: () => Promise<{ cancelled: boolean }>
  preview: (path: string, maxRows?: number) => Promise<Preview>
  onLog: (cb: (line: string) => void) => () => void
  onProgress: (cb: (p: { done: number; total: number }) => void) => () => void
  onDone: (cb: (result: RunResult) => void) => () => void
}

declare global {
  interface Window {
    fastag: FastagApi
  }
}
