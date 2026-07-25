// Where each CLI parameter appears in the UI.
//
// The manifest (params.generated.json) says what the parameters ARE; this says
// how they are PRESENTED. They are deliberately separate: the manifest is
// regenerated from the tool and must never be hand-edited, while this layout is
// a human judgement about which knobs a user reaches for.
//
// TOPP's own advanced="true" flag is not usable for this split -- it marks only
// TOPPBase boilerplate (log, debug, force, test, ...) and leaves every FASTag
// parameter, including deep internals like gap_penalty, marked non-advanced.

import manifest from './params.generated.json'

export interface ParamSpec {
  name: string
  type: string
  default: string | string[]
  description: string
  required: boolean
  toppAdvanced: boolean
  min?: number
  max?: number
  choices?: string[]
}

// UI overlays: the manifest is a faithful dump of the tool, but a few string
// params have a known small set of valid values the tool does not itself
// declare as restrictions. Turning them into a dropdown keeps the user from
// typing a rank the rollup can't resolve.
const CHOICES: Record<string, string[]> = {
  species_rank: ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']
}

export const PARAMS: ParamSpec[] = (manifest.params as ParamSpec[]).map((p) =>
  CHOICES[p.name] ? { ...p, choices: CHOICES[p.name] } : p
)
export const PARAM_BY_NAME = new Map(PARAMS.map((p) => [p.name, p]))

/// Shown in the main pane: the knobs that change a routine run. Species
/// detection lives here too -- enabling it is a first-class choice, not an
/// advanced tweak; only the taxonomy database paths are buried under Advanced.
export const CORE: string[] = [
  'tag_length',
  'fragment_tolerance',
  'fragment_tolerance_unit',
  'max_tags',
  'extension',
  'gaps',
  'deisotope',
  'max_evalue',
  'threads',
  'species',
  'species_rank',
  'species_min_len',
  'species_out'
]

/// Everything else worth exposing, grouped inside the Advanced accordion.
export const GROUPS: { title: string; params: string[] }[] = [
  { title: 'Output', params: ['proforma'] },
  { title: 'Peak selection', params: ['max_peaks', 'peaks_per_window'] },
  {
    title: 'Scoring & ranking',
    params: ['gap_penalty', 'orientation', 'isobaric_tolerance', 'min_filter_length']
  },
  { title: 'Modifications', params: ['fixed_modifications', 'variable_modifications'] },
  { title: 'Sequence database', params: ['fasta', 'out_spectra'] },
  { title: 'Subsampling', params: ['subsample_spectra', 'subsample_fraction', 'subsample_seed'] },
  // Only the taxonomy DB paths -- the enable/rank/output knobs are in CORE.
  { title: 'Taxonomy database', params: ['taxdb', 'taxonomy_nodes', 'taxonomy_names'] }
]

/// Deliberately not rendered, each for a stated reason. Listing them (rather
/// than defaulting to hide) is what lets the check below be exhaustive.
export const HIDDEN: Record<string, string> = {
  in: 'dedicated input picker',
  out: 'dedicated output field',
  progress: 'forced on; the GUI consumes the progress stream',
  version: 'TOPP boilerplate',
  log: 'TOPP boilerplate',
  debug: 'TOPP boilerplate',
  no_progress: 'TOPP boilerplate; unrelated to -progress',
  force: 'TOPP boilerplate',
  test: 'TOPP boilerplate'
}

/// The parameters actually shown, and so the only ones a run may send. HIDDEN
/// entries are not merely invisible: several are not settable at all (`-version`
/// records the tool's version in the INI and the CLI rejects it as an option),
/// so submitting the whole manifest aborts the run.
export const RENDERED: string[] = [...CORE, ...GROUPS.flatMap((g) => g.params)]

/// Every parameter must be placed somewhere. A new CLI option then shows up as
/// a loud failure instead of silently never appearing in the GUI -- the whole
/// point of generating the manifest from the tool.
export function unplacedParams(): string[] {
  const placed = new Set<string>([...CORE, ...GROUPS.flatMap((g) => g.params), ...Object.keys(HIDDEN)])
  return PARAMS.map((p) => p.name).filter((n) => !placed.has(n))
}
