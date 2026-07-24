import { contextBridge, ipcRenderer, IpcRendererEvent } from 'electron'

// The single, allowlisted bridge between renderer and main. The renderer gets
// exactly these functions and nothing else — no ipcRenderer, no Node.
const api = {
  probe: () => ipcRenderer.invoke('fastag:probe'),
  pickInput: () => ipcRenderer.invoke('dialog:pickInput'),
  pickOutput: (defaultPath?: string) => ipcRenderer.invoke('dialog:pickOutput', defaultPath),
  run: (params: unknown) => ipcRenderer.invoke('fastag:run', params),
  cancel: () => ipcRenderer.invoke('fastag:cancel'),
  preview: (path: string, maxRows?: number) => ipcRenderer.invoke('fastag:preview', path, maxRows),
  species: (path: string) => ipcRenderer.invoke('fastag:species', path),
  taxdbInfo: (explicit?: string) => ipcRenderer.invoke('fastag:taxdbInfo', explicit),
  loadSettings: () => ipcRenderer.invoke('settings:load'),
  saveLast: (values: unknown) => ipcRenderer.invoke('settings:saveLast', values),
  savePreset: (name: string, values: unknown) => ipcRenderer.invoke('settings:savePreset', name, values),
  deletePreset: (name: string) => ipcRenderer.invoke('settings:deletePreset', name),
  openTaxon: (taxid: number) => ipcRenderer.invoke('fastag:openTaxon', taxid),

  // Streamed events. Return an unsubscribe so the renderer can clean up.
  onLog: (cb: (line: string) => void) => {
    const h = (_e: IpcRendererEvent, line: string): void => cb(line)
    ipcRenderer.on('fastag:log', h)
    return () => ipcRenderer.removeListener('fastag:log', h)
  },
  onProgress: (cb: (p: { done: number; total: number }) => void) => {
    const h = (_e: IpcRendererEvent, p: { done: number; total: number }): void => cb(p)
    ipcRenderer.on('fastag:progress', h)
    return () => ipcRenderer.removeListener('fastag:progress', h)
  },
  onDone: (cb: (result: { ok: boolean; code: number | null; message?: string }) => void) => {
    const h = (_e: IpcRendererEvent, result: { ok: boolean; code: number | null; message?: string }): void =>
      cb(result)
    ipcRenderer.on('fastag:done', h)
    return () => ipcRenderer.removeListener('fastag:done', h)
  }
}

contextBridge.exposeInMainWorld('fastag', api)

export type FastagApi = typeof api
