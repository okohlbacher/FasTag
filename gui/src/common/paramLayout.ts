// Where each CLI parameter appears in the UI.
//
// The manifest (params.generated.json) says what the parameters ARE; this says
// how they are PRESENTED. They are deliberately separate: the manifest is
// regenerated from the tool and must never be hand-edited, while this layout is
// a human judgement about which knobs a user reaches for.
//
// TOPP's own advanced="true" flag is not usable for this split -- it marks only
// TOPPBase boilerplate (log, debug, force, test, ...) and leaves every FasTag
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

export const PARAMS: ParamSpec[] = manifest.params as ParamSpec[]
export const PARAM_BY_NAME = new Map(PARAMS.map((p) => [p.name, p]))

/// Shown in the main pane: the knobs that change a routine run.
export const CORE: string[] = [
  'tag_length',
  'fragment_tolerance',
  'fragment_tolerance_unit',
  'max_tags',
  'extension',
  'gaps',
  'deisotope',
  'max_evalue',
  'threads'
]

/// Everything else worth exposing, grouped inside the Advanced accordion.
export const GROUPS: { title: string; params: string[] }[] = [
  { title: 'Peak selection', params: ['max_peaks', 'peaks_per_window'] },
  {
    title: 'Scoring & ranking',
    params: ['gap_penalty', 'orientation', 'isobaric_tolerance', 'min_filter_length']
  },
  { title: 'Modifications', params: ['fixed_modifications', 'variable_modifications'] },
  { title: 'Sequence database', params: ['fasta', 'out_spectra'] },
  { title: 'Subsampling', params: ['subsample_spectra', 'subsample_fraction', 'subsample_seed'] },
  {
    title: 'Species detection',
    params: [
      'species',
      'taxdb',
      'taxonomy_nodes',
      'taxonomy_names',
      'species_out',
      'species_min_len',
      'species_rank'
    ]
  }
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
