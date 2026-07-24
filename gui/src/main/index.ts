import { app, shell, BrowserWindow, ipcMain, dialog } from 'electron'
import { join } from 'node:path'
import { probeBinary, runFastag, cancelRun, RunHandle, RunParams } from './fastag'
import { previewTsv } from './preview'

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
