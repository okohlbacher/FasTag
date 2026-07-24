import type { SpeciesReport, Taxon } from '../../preload/index.d'

// Genera that turn up in almost any LC-MS/MS run regardless of the sample:
// porcine trypsin, bovine serum, keratin from handling. Flagged so a user does
// not read a reagent as a finding. A hint, never a verdict -- these are also
// perfectly valid samples.
const CONTAMINANT: Record<string, string> = {
  Sus: 'porcine trypsin',
  Bos: 'bovine serum / BSA',
  Oryctolagus: 'common reagent source',
  Gallus: 'common reagent source',
  Ovis: 'common reagent source',
  Equus: 'common reagent source',
  Homo: 'keratin / handling (or the sample itself)'
}

function fmtQ(q: number): string {
  if (q === 0) return '<1e-300'
  if (q < 0.001) return q.toExponential(1)
  return q.toFixed(3)
}

interface Props {
  report: SpeciesReport | null
  enabled: boolean
  tooShort: { reach: number; k: number } | null
  contributed: number | null
  totalMs2: number | null
}

export default function SpeciesPanel({
  report,
  enabled,
  tooShort,
  contributed,
  totalMs2
}: Props): JSX.Element {
  if (!enabled) {
    return (
      <div className="empty">
        Species detection is off. Enable <code>species</code> under Advanced → Species detection.
      </div>
    )
  }

  if (tooShort) {
    return (
      <div className="warn big">
        <strong>No tag can reach the index.</strong>
        <p>
          Tags reach {tooShort.reach} residues but the index is keyed on {tooShort.k}-mers, so
          nothing can be looked up and the report would be empty. Raise <code>tag length</code> to{' '}
          {tooShort.k}, or use <code>extension</code> to reach it.
        </p>
      </div>
    )
  }

  if (!report) return <div className="empty">Species results appear here after a run.</div>
  if (report.empty) {
    return (
      <div className="empty">
        No taxon was supported by the tags. The run completed — this is a negative result, not an
        error.
      </div>
    )
  }

  const top = report.taxa[0]
  const maxObs = Math.max(...report.taxa.map((t) => t.observed), 1)
  const pctSpectra =
    contributed != null && totalMs2 ? Math.round((contributed / totalMs2) * 1000) / 10 : null

  return (
    <div className="species">
      <div className="headline">
        <div className="call">
          <span className="taxon">{top.name}</span>
          <span className="rank">{top.rank}</span>
        </div>
        <div className="sub">
          {top.observed.toLocaleString()} spectra · {top.enrichment}× enriched over background
        </div>
        <div className="caveat">
          {top.rank}-level call — a reduced reference set and short tags cannot resolve species.
        </div>
      </div>

      {contributed != null && (
        <div className="evidence">
          <strong>{contributed.toLocaleString()}</strong> spectra contributed taxon evidence
          {totalMs2 ? (
            <>
              {' '}
              of <strong>{totalMs2.toLocaleString()}</strong> MS2 ({pctSpectra}%)
            </>
          ) : null}
        </div>
      )}

      <div className="tablewrap">
        <table className="taxa">
          <thead>
            <tr>
              <th>taxon</th>
              <th className="num">observed</th>
              <th className="num">expected</th>
              <th className="num">enrich.</th>
              <th className="num">q</th>
            </tr>
          </thead>
          <tbody>
            {report.taxa.map((t: Taxon) => {
              const flag = CONTAMINANT[t.name]
              return (
                <tr key={t.taxid}>
                  <td>
                    <button
                      className="link"
                      title={`NCBI Taxonomy ${t.taxid}`}
                      onClick={() => window.fastag.openTaxon(t.taxid)}
                    >
                      {t.name}
                    </button>
                    {flag && (
                      <span className="flag" title={`Common contaminant: ${flag}`}>
                        ⚑
                      </span>
                    )}
                  </td>
                  <td className="num bar-cell">
                    <span className="minibar" style={{ width: `${(t.observed / maxObs) * 100}%` }} />
                    <span className="v">{t.observed.toLocaleString()}</span>
                  </td>
                  <td className="num dim">{t.expected.toLocaleString()}</td>
                  <td className="num">{t.enrichment}×</td>
                  <td className="num dim">{fmtQ(t.q)}</td>
                </tr>
              )
            })}
          </tbody>
        </table>
      </div>

      <p className="footnote">
        Ranked by significance, not by enrichment — the most enriched taxon is often not the right
        answer, because an over-represented proteome already carries a high background expectation.
        <br />
        q is a <em>ranking aid, not a calibrated FDR</em>: the per-k-mer background is a proxy and
        the chimeric-null problem applies. ⚑ marks genera common to most runs as reagents.
      </p>
    </div>
  )
}
