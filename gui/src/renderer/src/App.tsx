import { useEffect, useMemo, useRef, useState } from 'react'
import type { BinaryInfo, Preview, RunResult } from '../../preload/index.d'
import {
  CORE,
  GROUPS,
  PARAM_BY_NAME,
  RENDERED,
  unplacedParams,
  type ParamSpec
} from '../../common/paramLayout'
import ParamField, { type ParamValue } from './ParamField'
import SpeciesPanel from './SpeciesPanel'
import type { SpeciesReport } from '../../preload/index.d'
import logo from './assets/logo.svg'

// Where -species writes when -species_out is not given: <out> minus its final
// extension, plus .species.tsv. Mirrors the CLI so the GUI reads what it wrote.
function defaultSpeciesOut(out: string): string {
  const dot = out.lastIndexOf('.')
  return `${dot > 0 ? out.slice(0, dot) : out}.species.tsv`
}

interface Job {
  input: string
  out: string
  status: 'queued' | 'running' | 'done' | 'failed'
  tags?: number
}

function defaultOut(input: string): string {
  const dot = input.lastIndexOf('.')
  const stem = dot > 0 ? input.slice(0, dot) : input
  return `${stem}.tags.tsv`
}

// Initial form state IS the tool's own defaults, straight from the manifest.
function initialValues(): Record<string, ParamValue> {
  const v: Record<string, ParamValue> = {}
  for (const [name, spec] of PARAM_BY_NAME) {
    if (spec.type === 'string-list') v[name] = Array.isArray(spec.default) ? spec.default : []
    else if (spec.type === 'bool') v[name] = spec.default === 'true'
    else v[name] = String(spec.default ?? '')
  }
  return v
}

export default function App(): JSX.Element {
  const [bin, setBin] = useState<BinaryInfo | null>(null)
  const [input, setInput] = useState('')
  const [out, setOut] = useState('')
  const [values, setValues] = useState<Record<string, ParamValue>>(initialValues)
  const [advOpen, setAdvOpen] = useState(false)
  const [openGroups, setOpenGroups] = useState<Record<string, boolean>>({})
  const [running, setRunning] = useState(false)
  const [progress, setProgress] = useState<{ done: number; total: number } | null>(null)
  const [log, setLog] = useState<string[]>([])
  const [preview, setPreview] = useState<Preview | null>(null)
  const [species, setSpecies] = useState<SpeciesReport | null>(null)
  const [tab, setTab] = useState<'results' | 'species'>('results')
  const [taxdbK, setTaxdbK] = useState<number | null>(null)
  const [presets, setPresets] = useState<string[]>([])
  const [preset, setPreset] = useState('')
  const [savingName, setSavingName] = useState<string | null>(null)  // non-null = the 'save as' input is open
  const [jobs, setJobs] = useState<Job[]>([])
  const [batchRunning, setBatchRunning] = useState(false)
  // Awaiter for the run in flight, keyed by the main process's run id. Keying by
  // id (not a single slot) means a duplicated or late terminal event can only
  // ever resolve the run it belongs to -- never the next job's awaiter.
  const pending = useRef(new Map<number, (r: RunResult) => void>())
  // A terminal event can beat the invoke reply: a binary that vanishes after
  // probing makes the child emit 'error' before `run()` resolves with the run
  // id, so onDone would arrive before the resolver is registered. Park such a
  // result here; runOne claims it the moment it learns the id.
  const earlyDone = useRef(new Map<number, RunResult>())
  // Parsed from the CLI's own summary line, so the panel reports what the run
  // actually did rather than a number the GUI guessed.
  const [evidence, setEvidence] = useState<{ contributed: number | null; ms2: number | null }>({
    contributed: null,
    ms2: null
  })
  const logRef = useRef<HTMLPreElement>(null)

  // A parameter that exists in the CLI but nowhere in the layout would silently
  // be unreachable; surface it instead of hiding it.
  const unplaced = useMemo(() => unplacedParams(), [])

  const speciesOn = values['species'] === true
  // k comes from the index header, not an assumption: a differently-built index
  // would make a hardcoded 7 quietly wrong. The CLI warns too, but only after a
  // run has spent seconds and ~2 GB to produce an empty table.
  const reach = Number(values['tag_length'] || 0) + 2 * Number(values['extension'] || 0)
  const tooShort = speciesOn && taxdbK != null && reach < taxdbK ? { reach, k: taxdbK } : null

  useEffect(() => {
    window.fastag.probe().then(setBin)
    window.fastag.loadSettings().then((st) => {
      setPresets(Object.keys(st.presets).sort())
      // Restore the last session's parameters (but not in/out paths -- those are
      // per-run). Merge over defaults so a new CLI option added since last launch
      // still gets its default rather than undefined.
      if (st.lastUsed) setValues((v) => ({ ...v, ...(st.lastUsed as Record<string, ParamValue>) }))
    })
  }, [])

  // Re-read whenever the user points at a different index.
  const taxdbPath = String(values['taxdb'] || '')
  useEffect(() => {
    window.fastag.taxdbInfo(taxdbPath || undefined).then((i) => setTaxdbK(i ? i.k : null))
  }, [taxdbPath])

  useEffect(() => {
    const offLog = window.fastag.onLog((line) => {
      setLog((l) => [...l, line])
      // "Species: N spectra contributed taxon evidence; ..."
      const sp = /^Species:\s+(\d+) spectra contributed/.exec(line)
      if (sp) setEvidence((e) => ({ ...e, contributed: Number(sp[1]) }))
      // "42092 MS2 spectra, 1130228 tags"
      const ms2 = /^(\d+) MS2 spectra,/.exec(line)
      if (ms2) setEvidence((e) => ({ ...e, ms2: Number(ms2[1]) }))
    })
    const offProgress = window.fastag.onProgress((p) => setProgress(p))
    const offDone = window.fastag.onDone((r: RunResult) => {
      setProgress(null)
      setLog((l) => [...l, r.ok ? '— done —' : `— failed (${r.message ?? `exit ${r.code}`}) —`])
      // Hand the result to the awaiter for THIS run id. A terminal event with no
      // matching awaiter (a duplicate, or a run already settled by a failed
      // start) is dropped rather than resolving an unrelated job.
      const id = r.runId
      const resolve = id != null ? pending.current.get(id) : undefined
      if (id != null) pending.current.delete(id)
      if (resolve) resolve(r)
      else if (id != null) earlyDone.current.set(id, r) // beat the resolver; runOne will claim it
      if (pending.current.size === 0) setRunning(false)
    })
    return () => {
      offLog()
      offProgress()
      offDone()
    }
  }, [])

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight
  }, [log])

  const setValue = (name: string, v: ParamValue): void => setValues((s) => ({ ...s, [name]: v }))

  async function pickInput(): Promise<void> {
    const p = await window.fastag.pickInput()
    if (p) {
      setInput(p)
      setOut(defaultOut(p))
    }
  }

  async function addBatchFiles(): Promise<void> {
    const files = await window.fastag.pickInputs()
    if (!files.length) return
    setJobs((js) => {
      const have = new Set(js.map((j) => j.input))
      const add = files.filter((f) => !have.has(f)).map((f) => ({ input: f, out: defaultOut(f), status: 'queued' as const }))
      return [...js, ...add]
    })
  }

  async function pickFor(spec: ParamSpec): Promise<void> {
    const p =
      spec.type === 'output-file'
        ? await window.fastag.pickOutput(String(values[spec.name] || ''))
        : await window.fastag.pickInput()
    if (p) setValue(spec.name, p)
  }

  // Only the parameters the UI renders. The form seeds a value for every manifest
  // entry, but HIDDEN ones must never reach the command line: `-version` is
  // recorded in the INI yet rejected as an option, so sending the whole record
  // aborts the run with "Unknown option(s) '[-version]'".
  function submittedParams(): Record<string, ParamValue> {
    const s: Record<string, ParamValue> = {}
    for (const name of RENDERED) s[name] = values[name]
    return s
  }

  // Run one file and resolve when it finishes. Shared by single and batch runs.
  function runOne(inFile: string, outFile: string): Promise<RunResult> {
    return new Promise<RunResult>((resolve) => {
      setLog([])
      setProgress(null)
      setEvidence({ contributed: null, ms2: null }) // don't carry a prior run's counts
      setRunning(true)
      const settle = (r: RunResult): void => {
        if (pending.current.size === 0) setRunning(false)
        resolve(r)
      }
      window.fastag
        .run({ in: inFile, out: outFile, params: submittedParams() })
        .then((res) => {
          // Register the awaiter only once the main process hands back the run
          // id; the terminal event carries the same id (see onDone).
          if (res.started && res.runId != null) {
            const early = earlyDone.current.get(res.runId)
            if (early) {
              // The terminal event already fired before this reply landed.
              earlyDone.current.delete(res.runId)
              settle(early)
            } else {
              pending.current.set(res.runId, resolve)
            }
          } else {
            setLog([`could not start: ${res.reason ?? 'unknown'}`])
            settle({ ok: false, code: null, message: res.reason })
          }
        })
        .catch((err) => {
          // A rejected invoke would otherwise leave running=true forever.
          setLog([`could not start: ${err instanceof Error ? err.message : String(err)}`])
          settle({ ok: false, code: null, message: String(err) })
        })
    })
  }

  async function showResults(outFile: string): Promise<void> {
    setPreview(null)
    const rep0 = await window.fastag.preview(outFile, 200)
    setPreview(rep0)
    if (speciesOn) {
      const sp = await window.fastag.species(String(values['species_out'] || '') || defaultSpeciesOut(outFile))
      setSpecies(sp)
      if (sp && !sp.empty) setTab('species')
    }
  }

  async function run(): Promise<void> {
    if (!input || !out) return
    window.fastag.saveLast(values)
    const r = await runOne(input, out)
    if (r.ok) await showResults(out)
  }

  async function runBatch(): Promise<void> {
    if (jobs.length === 0 || batchRunning) return
    window.fastag.saveLast(values)
    setBatchRunning(true)
    // Snapshot the queue we're running and key status updates by input path (the
    // queue is dedup'd on input). A thrown preview/species must still release the
    // batch controls -- hence finally -- or Run stays disabled forever.
    const queue = jobs.filter((j) => j.status !== 'done')
    try {
      for (let i = 0; i < queue.length; i++) {
        const job = queue[i]
        setJobs((js) => js.map((j) => (j.input === job.input ? { ...j, status: 'running' } : j)))
        const r = await runOne(job.input, job.out)
        setJobs((js) => js.map((j) => (j.input === job.input ? { ...j, status: r.ok ? 'done' : 'failed' } : j)))
        if (r.ok && i === queue.length - 1) await showResults(job.out)  // preview the last
      }
    } finally {
      setBatchRunning(false)
    }
  }

  function resetDefaults(): void {
    setValues(initialValues())
  }

  async function applyPreset(name: string): Promise<void> {
    setPreset(name)
    if (!name) return
    const st = await window.fastag.loadSettings()
    const pv = st.presets[name]
    if (pv) setValues((v) => ({ ...v, ...(pv as Record<string, ParamValue>) }))
  }
  async function commitPreset(): Promise<void> {
    const name = (savingName ?? '').trim()
    if (!name) { setSavingName(null); return }
    await window.fastag.savePreset(name, values)   // window.prompt is unsupported in Electron; inline input instead
    const st = await window.fastag.loadSettings()
    setPresets(Object.keys(st.presets).sort())
    setPreset(name)
    setSavingName(null)
  }
  async function removePreset(): Promise<void> {
    if (!preset) return
    await window.fastag.deletePreset(preset)
    const st = await window.fastag.loadSettings()
    setPresets(Object.keys(st.presets).sort())
    setPreset('')
  }

  const canRun = bin?.ok && !!input && !!out && !running && !batchRunning
  const field = (name: string): JSX.Element | null => {
    const spec = PARAM_BY_NAME.get(name)
    if (!spec) return null
    return (
      <ParamField
        key={name}
        spec={spec}
        value={values[name]}
        onChange={(v) => setValue(name, v)}
        onPickFile={pickFor}
      />
    )
  }

  return (
    <div className="app">
      <header>
        <img className="logo" src={logo} alt="" width={26} height={26} />
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
            <div className="row tight">
              <button className="secondary" onClick={pickInput}>
                Choose file…
              </button>
              <button className="secondary" onClick={addBatchFiles} disabled={batchRunning || running} title="Queue several files to run in turn">
                Add batch…
              </button>
            </div>
            {input && jobs.length === 0 && <div className="path">{input}</div>}
          </div>

          {jobs.length > 0 && (
            <div className="field batch">
              <label>
                Batch queue
                <span className="count">
                  {jobs.filter((j) => j.status === 'done').length}/{jobs.length} done
                </span>
              </label>
              <ul className="jobs">
                {jobs.map((j, i) => (
                  <li key={j.input} className={`job ${j.status}`}>
                    <span className="st" aria-hidden>
                      {j.status === 'done' ? '✓' : j.status === 'failed' ? '✕' : j.status === 'running' ? '⟳' : '·'}
                    </span>
                    <span className="jn" title={j.input}>{j.input.split('/').pop()}</span>
                    {!batchRunning && (
                      <button className="rm" title="remove" onClick={() => setJobs((js) => js.filter((_, k) => k !== i))}>
                        ×
                      </button>
                    )}
                  </li>
                ))}
              </ul>
              <div className="row tight">
                <button onClick={runBatch} disabled={batchRunning || !bin?.ok}>
                  {batchRunning ? 'Running batch…' : `Run batch (${jobs.filter((j) => j.status !== 'done').length})`}
                </button>
                <button className="secondary slim" onClick={() => setJobs([])} disabled={batchRunning}>
                  Clear
                </button>
              </div>
            </div>
          )}

          {jobs.length === 0 && (
            <div className="field">
              <label htmlFor="out">Output tags (TSV)</label>
              <input id="out" value={out} onChange={(e) => setOut(e.target.value)} placeholder="tags.tsv" />
            </div>
          )}

          <div className="sep" />

          {CORE.map(field)}

          <div className="presets">
            {savingName === null ? (
              <>
                <select value={preset} onChange={(e) => applyPreset(e.target.value)} title="Load a saved preset">
                  <option value="">Presets…</option>
                  {presets.map((n) => (
                    <option key={n} value={n}>{n}</option>
                  ))}
                </select>
                <button className="secondary slim" onClick={() => setSavingName(preset || '')}>Save as…</button>
                <button className="secondary slim" onClick={removePreset} disabled={!preset}>Delete</button>
              </>
            ) : (
              <>
                <input
                  autoFocus
                  value={savingName}
                  placeholder="preset name"
                  onChange={(e) => setSavingName(e.target.value)}
                  onKeyDown={(e) => { if (e.key === 'Enter') commitPreset(); if (e.key === 'Escape') setSavingName(null) }}
                />
                <button className="slim" onClick={commitPreset}>Save</button>
                <button className="secondary slim" onClick={() => setSavingName(null)}>Cancel</button>
              </>
            )}
          </div>

          <div className="adv">
            <button className="disclosure" onClick={() => setAdvOpen((o) => !o)} aria-expanded={advOpen}>
              <span className={`chev ${advOpen ? 'open' : ''}`}>▸</span> Advanced
              <span className="count">
                {GROUPS.reduce((n, g) => n + g.params.length, 0)} settings
              </span>
            </button>
            {advOpen && (
              <div className="adv-body">
                {GROUPS.map((g) => {
                  const open = openGroups[g.title] ?? false
                  return (
                    <div className="group" key={g.title}>
                      <button
                        className="disclosure sub"
                        onClick={() => setOpenGroups((s) => ({ ...s, [g.title]: !open }))}
                        aria-expanded={open}
                      >
                        <span className={`chev ${open ? 'open' : ''}`}>▸</span> {g.title}
                        <span className="count">{g.params.length}</span>
                      </button>
                      {open && <div className="group-body">{g.params.map(field)}</div>}
                    </div>
                  )
                })}
                {unplaced.length > 0 && (
                  <p className="warn">
                    {unplaced.length} CLI parameter(s) missing from the UI layout: {unplaced.join(', ')}
                  </p>
                )}
                <button className="secondary slim" onClick={resetDefaults}>
                  Reset all to CLI defaults
                </button>
              </div>
            )}
          </div>

          <div className="row actions">
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
                    {Math.round((progress.done / progress.total) * 100)}% ·{' '}
                    {progress.done.toLocaleString()}/{progress.total.toLocaleString()} spectra
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

          <div className="tabs" role="tablist">
            <button
              role="tab"
              aria-selected={tab === 'results'}
              className={tab === 'results' ? 'on' : ''}
              onClick={() => setTab('results')}
            >
              Tags
              {preview && <span className="count">{preview.shown}{preview.truncated ? '+' : ''}</span>}
            </button>
            <button
              role="tab"
              aria-selected={tab === 'species'}
              className={tab === 'species' ? 'on' : ''}
              onClick={() => setTab('species')}
            >
              Species
              {species && !species.empty && <span className="count">{species.taxa.length}</span>}
            </button>
          </div>

          {tab === 'species' ? (
            <div className="tablewrap">
              <SpeciesPanel
                report={species}
                enabled={speciesOn}
                tooShort={tooShort}
                contributed={evidence.contributed}
                totalMs2={evidence.ms2}
              />
            </div>
          ) : (
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
          )}
        </section>
      </main>
    </div>
  )
}
