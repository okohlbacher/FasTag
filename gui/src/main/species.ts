// Read the ranked-taxa TSV that -species writes.
//
// Unlike the tag table this file is small (one row per taxon), so it is read
// whole. Columns, from FasTag: rank, taxid, name, observed, expected,
// enrichment, log_pvalue, qvalue.

import { existsSync, readFileSync, openSync, readSync, closeSync } from 'node:fs'
import { join, dirname } from 'node:path'

export interface TaxdbInfo {
  path: string
  k: number
  kmers: number
}

// Read k straight out of the index header rather than assuming it.
//
// Layout (TaxIndex::save): "FTXI" | uint32 version | uint32 k | uint64 nKmers.
// The GUI needs k to tell the user, BEFORE a run, that their tag_length cannot
// reach the index -- the failure mode otherwise is a successful-looking empty
// report after a ~2 GB index load. Hardcoding 7 would silently mislead the
// moment someone builds an index at another k.
export function readTaxdbInfo(path: string): TaxdbInfo | null {
  if (!path || !existsSync(path)) return null
  let fd: number | null = null
  try {
    fd = openSync(path, 'r')
    const buf = Buffer.alloc(20)
    if (readSync(fd, buf, 0, 20, 0) < 20) return null
    const magic = buf.toString('latin1', 0, 4)
    // Both versions store k at offset 8. The k-mer count moved: v1 (FTXI) put n
    // at 12; v2 (FTX2) put n_taxa at 12 and n_kmers at 16.
    if (magic !== 'FTXI' && magic !== 'FTX2') return null
    const kmers = magic === 'FTX2' ? Number(buf.readBigUInt64LE(16)) : Number(buf.readBigUInt64LE(12))
    return { path, k: buf.readUInt32LE(8), kmers }
  } catch {
    return null
  } finally {
    if (fd !== null) closeSync(fd)
  }
}

// The same search order the CLI uses (see taxonomyDir_ in FasTag.cpp), so the
// GUI reports on the index the run will actually load.
export function bundledTaxdb(binaryPath: string): string | null {
  const env = process.env.FASTAG_TAXONOMY_DIR
  const dirs: string[] = []
  if (env) dirs.push(env)
  const bin = dirname(binaryPath)
  dirs.push(join(bin, '..', 'share', 'FasTag', 'taxonomy'))
  dirs.push(join(bin, 'share-FasTag-taxonomy'))
  dirs.push(join(bin, 'taxonomy'))
  for (const d of dirs) {
    const p = join(d, 'tax_k7.taxdb')
    if (existsSync(p)) return p
  }
  return null
}

export interface Taxon {
  rank: string
  taxid: number
  name: string
  observed: number
  expected: number
  enrichment: number
  logP: number
  q: number
}

export interface SpeciesReport {
  path: string
  taxa: Taxon[]
  /// True when the file exists but holds only a header -- the run asked for
  /// species detection and nothing supported a taxon. Distinct from "no file".
  empty: boolean
}

// FasTag writes enrichment as e.g. "1.4x"; keep the number.
function num(s: string): number {
  const v = parseFloat(String(s).replace(/x$/i, ''))
  return Number.isFinite(v) ? v : 0
}

export function readSpecies(path: string): SpeciesReport | null {
  if (!path || !existsSync(path)) return null
  const lines = readFileSync(path, 'utf8').split('\n').filter((l) => l.trim() !== '')
  if (lines.length === 0) return null

  const head = lines[0].split('\t')
  const col = (n: string): number => head.indexOf(n)
  const iRank = col('rank')
  const iTaxid = col('taxid')
  const iName = col('name')
  const iObs = col('observed')
  const iExp = col('expected')
  const iEnr = col('enrichment')
  const iLogP = col('log_pvalue')
  const iQ = col('qvalue')

  const taxa: Taxon[] = []
  for (const line of lines.slice(1)) {
    const f = line.split('\t')
    taxa.push({
      rank: f[iRank] ?? '',
      taxid: Number(f[iTaxid] ?? 0),
      name: f[iName] ?? '',
      observed: num(f[iObs] ?? '0'),
      expected: num(f[iExp] ?? '0'),
      enrichment: num(f[iEnr] ?? '0'),
      logP: num(f[iLogP] ?? '0'),
      q: num(f[iQ] ?? '1')
    })
  }

  // Rank by significance, NOT by enrichment. On a human sample Homo is the true
  // answer at 1.4x while a bit-part taxon shows 4.0x -- human dominates the
  // index, so its background expectation is already high. Sorting by enrichment
  // confidently returns the wrong organism.
  taxa.sort((a, b) => a.logP - b.logP || b.observed - a.observed)

  return { path, taxa, empty: taxa.length === 0 }
}
