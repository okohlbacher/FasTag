// Extract the FasTag parameter contract from the tool's own -write_ini output.
//
// The CLI is the source of truth for names, types, defaults and restrictions.
// Hand-copying them into the UI is how a GUI drifts from its CLI -- the first
// draft of this app hardcoded fragment_tolerance=0.02 with no unit, which meant
// 0.02 ppm and produced zero tags on every file. This script makes that class of
// mistake impossible to introduce by hand.
//
//   node scripts/gen-params.mjs [path/to/FasTag] > src/common/params.generated.json
//
// The binary defaults to the one the app itself resolves (resources/fastag), so
// running it under dev.sh picks up the right env.

import { execFileSync } from 'node:child_process'
import { mkdtempSync, readFileSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'

const bin = process.argv[2] || process.env.FASTAG_BIN || 'resources/fastag/bin/FasTag'

const dir = mkdtempSync(join(tmpdir(), 'fastag-ini-'))
const ini = join(dir, 'fastag.ini')
try {
  execFileSync(bin, ['-write_ini', ini], { stdio: ['ignore', 'ignore', 'pipe'] })
} catch (e) {
  console.error(`failed to run ${bin} -write_ini: ${e.message}`)
  process.exit(1)
}
const xml = readFileSync(ini, 'utf8')
rmSync(dir, { recursive: true, force: true })

const unescape = (s) =>
  s
    .replace(/&quot;/g, '"')
    .replace(/&apos;/g, "'")
    .replace(/&lt;/g, '<')
    .replace(/&gt;/g, '>')
    .replace(/&amp;/g, '&')

const attrs = (tag) => {
  const out = {}
  for (const m of tag.matchAll(/(\w+)="([^"]*)"/g)) out[m[1]] = unescape(m[2])
  return out
}

// restrictions carry two different meanings depending on type:
//   numeric -> "min:max", either side may be empty ("1:" = min 1, no max)
//   string  -> comma-separated allowed values
const parseRestrictions = (type, r) => {
  if (!r) return {}
  if (type === 'int' || type === 'double') {
    const [lo, hi] = r.split(':')
    const res = {}
    if (lo !== '' && lo !== undefined) res.min = Number(lo)
    if (hi !== '' && hi !== undefined) res.max = Number(hi)
    return res
  }
  return { choices: r.split(',') }
}

const params = []

for (const m of xml.matchAll(/<ITEM\s([^>]*?)\/>/g)) {
  const a = attrs(m[1])
  if (!a.name) continue
  params.push({
    name: a.name,
    type: a.type,
    default: a.value ?? '',
    description: a.description ?? '',
    required: a.required === 'true',
    toppAdvanced: a.advanced === 'true',
    ...parseRestrictions(a.type, a.restrictions)
  })
}

// <ITEMLIST name=... > <LISTITEM value="..."/> ... </ITEMLIST>
for (const m of xml.matchAll(/<ITEMLIST\s([^>]*?)>([\s\S]*?)<\/ITEMLIST>/g)) {
  const a = attrs(m[1])
  if (!a.name) continue
  const items = [...m[2].matchAll(/<LISTITEM\s+value="([^"]*)"/g)].map((x) => unescape(x[1]))
  params.push({
    name: a.name,
    type: 'string-list',
    default: items,
    description: a.description ?? '',
    required: a.required === 'true',
    toppAdvanced: a.advanced === 'true'
  })
}

params.sort((x, y) => x.name.localeCompare(y.name))
process.stdout.write(JSON.stringify({ params }, null, 2) + '\n')
