import { app, shell, BrowserWindow, ipcMain, dialog } from 'electron'
import { join } from 'node:path'
import { probeBinary, resolveBinary, runFastag, cancelRun, RunHandle, RunParams } from './fastag'
import { previewTsv } from './preview'
import { readSpecies, readTaxdbInfo, bundledTaxdb } from './species'
import { loadSettings, saveLastUsed, savePreset, deletePreset } from './settings'

// One in-flight run for the P0/P1 skeleton (batch queue comes later).
let currentRun: RunHandle | null = null

function createWindow(): void {
  const win = new BrowserWindow({
    width: 1100,
    height: 760,
    show: false,
    title: 'FasTag',
    webPreferences: {
      preload: join(__dirname, '../preload/index.js'),
      // Hardening: renderer gets no Node, no remote, isolated context.
      sandbox: true,
      contextIsolation: true,
      nodeIntegration: false
    }
  })

  win.on('ready-to-show', () => win.show())

  // External links open in the OS browser, never in-app.
  win.webContents.setWindowOpenHandler(({ url }) => {
    shell.openExternal(url)
    return { action: 'deny' }
  })

  if (process.env['ELECTRON_RENDERER_URL']) {
    win.loadURL(process.env['ELECTRON_RENDERER_URL'])
  } else {
    win.loadFile(join(__dirname, '../renderer/index.html'))
  }
}

// ---- IPC: every handler is explicit and allowlisted -----------------------

ipcMain.handle('fastag:probe', () => probeBinary())

ipcMain.handle('dialog:pickInput', async () => {
  const r = await dialog.showOpenDialog({
    title: 'Select spectra file',
    properties: ['openFile'],
    filters: [
      { name: 'Spectra', extensions: ['mzML', 'mzml', 'mzpeak', 'mzPeak'] },
      { name: 'All files', extensions: ['*'] }
    ]
  })
  return r.canceled ? null : r.filePaths[0]
})

ipcMain.handle('dialog:pickOutput', async (_e, defaultPath?: string) => {
  const r = await dialog.showSaveDialog({
    title: 'Save tags as',
    defaultPath,
    filters: [{ name: 'Tags (TSV)', extensions: ['tsv'] }]
  })
  return r.canceled ? null : r.filePath
})

ipcMain.handle('fastag:run', (e, params: RunParams) => {
  if (currentRun) return { started: false, reason: 'a run is already in progress' }
  const wc = e.sender
  // The GUI always wants machine-readable progress; force the flag on.
  currentRun = runFastag(
    { ...params, progress: true },
    {
      onLog: (line) => wc.send('fastag:log', line),
      onProgress: (done, total) => wc.send('fastag:progress', { done, total }),
      onError: (message) => wc.send('fastag:done', { ok: false, code: null, message }),
      onExit: (code) => {
        currentRun = null
        wc.send('fastag:done', { ok: code === 0, code })
      }
    }
  )
  return { started: true }
})

ipcMain.handle('fastag:cancel', () => {
  if (currentRun) {
    cancelRun(currentRun)
    return { cancelled: true }
  }
  return { cancelled: false }
})

ipcMain.handle('fastag:preview', (_e, path: string, maxRows?: number) => previewTsv(path, maxRows))

ipcMain.handle('settings:load', () => loadSettings())
ipcMain.handle('settings:saveLast', (_e, values: Record<string, unknown>) => saveLastUsed(values))
ipcMain.handle('settings:savePreset', (_e, name: string, values: Record<string, unknown>) => savePreset(name, values))
ipcMain.handle('settings:deletePreset', (_e, name: string) => deletePreset(name))

ipcMain.handle('fastag:species', (_e, path: string) => readSpecies(path))

// k of the index a run would load, so the UI can warn before the run rather
// than after a ~2 GB load produces an empty report.
ipcMain.handle('fastag:taxdbInfo', (_e, explicit?: string) => {
  const p = explicit && explicit.trim() !== '' ? explicit : bundledTaxdb(resolveBinary().bin)
  return p ? readTaxdbInfo(p) : null
})

// Open the NCBI Taxonomy page for a taxid. Built here from a fixed template --
// the renderer passes a number, never a URL, so a crafted value cannot redirect
// the user somewhere else.
ipcMain.handle('fastag:openTaxon', (_e, taxid: number) => {
  const id = Number(taxid)
  if (!Number.isInteger(id) || id <= 0) return false
  shell.openExternal(`https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=${id}`)
  return true
})

// ---- lifecycle ------------------------------------------------------------

app.whenReady().then(() => {
  createWindow()
  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow()
  })
})

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit()
})
