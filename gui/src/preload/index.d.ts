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
}

export interface FastagApi {
  probe: () => Promise<BinaryInfo>
  pickInput: () => Promise<string | null>
  pickOutput: (defaultPath?: string) => Promise<string | null>
  run: (params: Record<string, unknown>) => Promise<{ started: boolean; reason?: string }>
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
