// Bounded TSV preview for the results pane (P0/P1).
//
// A full run is 100k-1M+ rows; the P0 preview reads only the header plus the
// first N data rows by streaming and stopping early, so a huge file never
// loads into memory. The million-row browser (DuckDB in a utility process)
// replaces this in a later phase — this is the honest capped stand-in.

import { createReadStream } from 'node:fs'
import { createInterface } from 'node:readline'

export interface Preview {
  header: string[]
  rows: string[][]
  truncated: boolean
  shown: number
}

export async function previewTsv(path: string, maxRows = 200): Promise<Preview> {
  return new Promise((resolve, reject) => {
    const rl = createInterface({ input: createReadStream(path, { encoding: 'utf8' }), crlfDelay: Infinity })
    let header: string[] | null = null
    const rows: string[][] = []
    let truncated = false

    rl.on('line', (line) => {
      if (header === null) {
        header = line.split('\t')
        return
      }
      if (rows.length < maxRows) {
        rows.push(line.split('\t'))
      } else {
        truncated = true
        rl.close()
      }
    })
    rl.on('close', () => resolve({ header: header ?? [], rows, truncated, shown: rows.length }))
    rl.on('error', reject)
  })
}
