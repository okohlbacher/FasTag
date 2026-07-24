import { app, shell, BrowserWindow, ipcMain, dialog } from 'electron'
import { join } from 'node:path'
import { probeBinary, resolveBinary, runFastag, cancelRun, RunHandle, RunParams } from './fastag'
import { previewTsv } from './preview'
import { readSpecies, readTaxdbInfo, bundledTaxdb } from './species'
import { loadSettings, saveLastUsed, savePreset, deletePreset } from './settings'

// One in-flight run (the batch queue runs sequentially, one process at a time).
let currentRun: RunHandle | null = null
// Monotonic id so a terminal event is correlated to the run that produced it:
// the renderer only resolves the awaiter whose id matches, never a stray one.
let runSeq = 0

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

  // Killing the window (close or reload) must not orphan a running CLI process.
  const killOrphan = (): void => {
    if (currentRun) {
      cancelRun(currentRun)
      currentRun = null
    }
  }
  win.on('closed', killOrphan)
  win.webContents.on('did-start-navigation', (_e, _url, isInPlace, isMainFrame) => {
    if (isMainFrame && !isInPlace) killOrphan() // a reload/navigation, not a hash change
  })

  // External links open in the OS browser, never in-app -- and only real web
  // links: never hand file:, javascript:, or other schemes to the OS opener.
  win.webContents.setWindowOpenHandler(({ url }) => {
    if (/^https?:\/\//i.test(url)) shell.openExternal(url)
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

ipcMain.handle('dialog:pickInputs', async () => {
  const r = await dialog.showOpenDialog({
    title: 'Select spectra files for batch',
    properties: ['openFile', 'multiSelections'],
    filters: [
      { name: 'Spectra', extensions: ['mzML', 'mzml', 'mzpeak', 'mzPeak'] },
      { name: 'All files', extensions: ['*'] }
    ]
  })
  return r.canceled ? [] : r.filePaths
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
  const runId = ++runSeq
  // Sending to a destroyed WebContents throws "Object has been destroyed" and
  // would take down the main process; a run whose window closed/reloaded still
  // emits log + terminal callbacks, so every send is guarded.
  const send = (channel: string, payload: unknown): void => {
    if (!wc.isDestroyed()) wc.send(channel, payload)
  }
  // The GUI always wants machine-readable progress; force the flag on.
  let handle: RunHandle | null = null
  // Clear currentRun only if it is still OURS: a stale run that ignored SIGTERM
  // can exit after a reload started a new run, and an unconditional null would
  // steal the new run's ownership (cancel could no longer stop it).
  const clearIfMine = (): void => {
    if (currentRun === handle) currentRun = null
  }
  handle = runFastag(
    { ...params, progress: true },
    {
      onLog: (line) => send('fastag:log', line),
      onProgress: (done, total) => send('fastag:progress', { done, total }),
      onError: (message) => {
        clearIfMine()
        send('fastag:done', { runId, ok: false, code: null, message })
      },
      onExit: (code) => {
        clearIfMine()
        send('fastag:done', { runId, ok: code === 0, code })
      }
    }
  )
  currentRun = handle
  return { started: true, runId }
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

// A second launch should focus the existing window, not spin up a rival main
// process that fights over the same settings file and run state.
if (!app.requestSingleInstanceLock()) {
  app.quit()
} else {
  app.on('second-instance', () => {
    const win = BrowserWindow.getAllWindows()[0]
    if (win) {
      if (win.isMinimized()) win.restore()
      win.focus()
    } else {
      // macOS keeps the app alive with no window; a relaunch should show one.
      createWindow()
    }
  })
  app.whenReady().then(() => {
    createWindow()
    app.on('activate', () => {
      if (BrowserWindow.getAllWindows().length === 0) createWindow()
    })
  })
}

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit()
})
