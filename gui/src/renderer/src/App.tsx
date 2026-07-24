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

  // Latest run context for the IPC subscription below, which registers once.
  const runCtx = useRef({ out: '', speciesOn: false, speciesOut: '' })
  runCtx.current = { out, speciesOn, speciesOut: String(values['species_out'] || '') }

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
      setRunning(false)
      setProgress(null)
      setLog((l) => [...l, r.ok ? '— done —' : `— failed (${r.message ?? `exit ${r.code}`}) —`])
      const ctx = runCtx.current
      if (r.ok && ctx.out) {
        window.fastag.preview(ctx.out, 200).then(setPreview)
        if (ctx.speciesOn) {
          window.fastag.species(ctx.speciesOut || defaultSpeciesOut(ctx.out)).then((rep) => {
            setSpecies(rep)
            if (rep && !rep.empty) setTab('species')
          })
        }
      }
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

  async function pickFor(spec: ParamSpec): Promise<void> {
    const p =
      spec.type === 'output-file'
        ? await window.fastag.pickOutput(String(values[spec.name] || ''))
        : await window.fastag.pickInput()
    if (p) setValue(spec.name, p)
  }

  async function run(): Promise<void> {
    if (!input || !out) return
    setLog([])
    setPreview(null)
    setProgress(null)
    setRunning(true)
    window.fastag.saveLast(values)   // remember this run's parameters for next launch
    // Only the parameters the UI actually renders. The form seeds a value for
    // every manifest entry, but HIDDEN ones must never reach the command line:
    // `-version` is recorded in the INI yet rejected as an option, so sending
    // the whole record aborts the run with "Unknown option(s) '[-version]'".
    const submitted: Record<string, ParamValue> = {}
    for (const name of RENDERED) submitted[name] = values[name]
    const res = await window.fastag.run({ in: input, out, params: submitted })
    if (!res.started) {
      setRunning(false)
      setLog([`could not start: ${res.reason ?? 'unknown'}`])
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

  const canRun = bin?.ok && !!input && !!out && !running
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
            <button className="secondary" onClick={pickInput}>
              Choose file…
            </button>
            {input && <div className="path">{input}</div>}
          </div>

          <div className="field">
            <label htmlFor="out">Output tags (TSV)</label>
            <input id="out" value={out} onChange={(e) => setOut(e.target.value)} placeholder="tags.tsv" />
          </div>

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
