import { useEffect, useRef, useState } from 'react'
import type { BinaryInfo, Preview, RunResult } from '../../preload/index.d'

// Derive a default output path next to the input: foo.mzML -> foo.tags.tsv
function defaultOut(input: string): string {
  const dot = input.lastIndexOf('.')
  const stem = dot > 0 ? input.slice(0, dot) : input
  return `${stem}.tags.tsv`
}

export default function App(): JSX.Element {
  const [bin, setBin] = useState<BinaryInfo | null>(null)
  const [input, setInput] = useState('')
  const [out, setOut] = useState('')
  // Defaults mirror the FasTag CLI's own defaults (from -write_ini): the GUI
  // must not invent its own, or a run with "defaults" silently differs from the
  // CLI. tag_length=3, fragment_tolerance=20 ppm, max_tags=50.
  const [tagLength, setTagLength] = useState(3)
  const [tol, setTol] = useState(20)
  const [tolUnit, setTolUnit] = useState<'ppm' | 'Da'>('ppm')
  const [maxTags, setMaxTags] = useState(50)
  const [running, setRunning] = useState(false)
  const [progress, setProgress] = useState<{ done: number; total: number } | null>(null)
  const [log, setLog] = useState<string[]>([])
  const [preview, setPreview] = useState<Preview | null>(null)
  const logRef = useRef<HTMLPreElement>(null)

  // Probe the binary once on mount — the architecture canary.
  useEffect(() => {
    window.fastag.probe().then(setBin)
  }, [])

  // Subscribe to streamed log + completion for the lifetime of the component.
  useEffect(() => {
    const offLog = window.fastag.onLog((line) => setLog((l) => [...l, line]))
    const offProgress = window.fastag.onProgress((p) => setProgress(p))
    const offDone = window.fastag.onDone((r: RunResult) => {
      setRunning(false)
      setProgress(null)
      setLog((l) => [...l, r.ok ? '— done —' : `— failed (${r.message ?? `exit ${r.code}`}) —`])
      if (r.ok && out) window.fastag.preview(out, 200).then(setPreview)
    })
    return () => {
      offLog()
      offProgress()
      offDone()
    }
  }, [out])

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight
  }, [log])

  async function pickInput(): Promise<void> {
    const p = await window.fastag.pickInput()
    if (p) {
      setInput(p)
      setOut(defaultOut(p))
    }
  }

  async function run(): Promise<void> {
    if (!input || !out) return
    setLog([])
    setPreview(null)
    setProgress(null)
    setRunning(true)
    const res = await window.fastag.run({
      in: input,
      out,
      tag_length: tagLength,
      fragment_tolerance: tol,
      fragment_tolerance_unit: tolUnit,
      max_tags: maxTags
    })
    if (!res.started) {
      setRunning(false)
      setLog([`could not start: ${res.reason ?? 'unknown'}`])
    }
  }

  const canRun = bin?.ok && !!input && !!out && !running

  return (
    <div className="app">
      <header>
        <h1>FasTag</h1>
        {bin && (
          <span className={`badge ${bin.ok ? 'ok' : 'bad'}`} title={bin.bin}>
            {bin.ok ? bin.detail : `binary not runnable: ${bin.detail}`}
          </span>
        )}
      </header>

      <main>
        <section className="controls">
          <div className="field">
            <label>Input spectra (mzML / mzPeak)</label>
            <button className="secondary" onClick={pickInput}>
              Choose file…
            </button>
            {input && <div className="path">{input}</div>}
          </div>

          <div className="field">
            <label>Output tags (TSV)</label>
            <input value={out} onChange={(e) => setOut(e.target.value)} placeholder="tags.tsv" />
          </div>

          <div className="row">
            <div className="field">
              <label>Tag length</label>
              <input
                type="number"
                min={2}
                max={10}
                value={tagLength}
                onChange={(e) => setTagLength(Number(e.target.value))}
              />
            </div>
            <div className="field">
              <label>Max tags / spectrum</label>
              <input
                type="number"
                min={1}
                max={100}
                value={maxTags}
                onChange={(e) => setMaxTags(Number(e.target.value))}
              />
            </div>
          </div>

          <div className="field">
            <label>Fragment tolerance</label>
            <div className="row">
              <input
                type="number"
                step={0.001}
                min={0}
                value={tol}
                onChange={(e) => setTol(Number(e.target.value))}
              />
              <select value={tolUnit} onChange={(e) => setTolUnit(e.target.value as 'ppm' | 'Da')}>
                <option value="ppm">ppm</option>
                <option value="Da">Da</option>
              </select>
            </div>
          </div>

          <div className="row">
            <button onClick={run} disabled={!canRun}>
              {running ? 'Running…' : 'Run'}
            </button>
            <button className="secondary" onClick={() => window.fastag.cancel()} disabled={!running}>
              Cancel
            </button>
          </div>
        </section>

        <section className="results">
          {running && (
            <div className="progress">
              {progress && progress.total > 0 ? (
                <>
                  <div className="bar">
                    <div
                      className="fill"
                      style={{ width: `${Math.round((progress.done / progress.total) * 100)}%` }}
                    />
                  </div>
                  <span className="pct">
                    {Math.round((progress.done / progress.total) * 100)}% · {progress.done.toLocaleString()}/
                    {progress.total.toLocaleString()} spectra
                  </span>
                </>
              ) : (
                <>
                  <div className="bar indeterminate">
                    <div className="fill" />
                  </div>
                  <span className="pct">
                    {progress ? `${progress.done.toLocaleString()} spectra` : 'starting…'}
                  </span>
                </>
              )}
            </div>
          )}
          <pre className="log" ref={logRef}>
            {log.length ? log.join('\n') : 'Log output appears here.'}
          </pre>
          <div className="tablewrap">
            {preview ? (
              preview.rows.length ? (
                <>
                  <table>
                    <thead>
                      <tr>
                        {preview.header.map((h, i) => (
                          <th key={i}>{h}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {preview.rows.map((r, i) => (
                        <tr key={i}>
                          {r.map((c, j) => (
                            <td key={j}>{c}</td>
                          ))}
                        </tr>
                      ))}
                    </tbody>
                  </table>
                  {preview.truncated && (
                    <div className="trunc">
                      Showing first {preview.shown} rows (preview cap — full browser lands in a later phase).
                    </div>
                  )}
                </>
              ) : (
                <div className="empty">Run produced no tags.</div>
              )
            ) : (
              <div className="empty">Results preview appears here after a run.</div>
            )}
          </div>
        </section>
      </main>
    </div>
  )
}
