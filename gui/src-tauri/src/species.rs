// Read the ranked-taxa TSV that -species writes, and read k from the index
// header so the UI can warn before a run that tag_length can't reach the index.
// Ported from gui/src/main/species.ts (and fixes a latent 20-vs-24 byte read
// that made the Electron build return null for every FTX2 index).

use serde::Serialize;
use std::io::Read;
use std::path::{Path, PathBuf};

use tauri::AppHandle;

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Taxon {
    rank: String,
    taxid: i64,
    name: String,
    observed: f64,
    expected: f64,
    enrichment: f64,
    log_p: f64,
    q: f64,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SpeciesReport {
    path: String,
    taxa: Vec<Taxon>,
    empty: bool,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct TaxdbInfo {
    path: String,
    k: u32,
    kmers: u64,
}

// FASTag writes enrichment as e.g. "1.4x"; keep the number.
fn num(s: &str) -> f64 {
    let t = s.trim().trim_end_matches(['x', 'X']);
    t.parse::<f64>().ok().filter(|v| v.is_finite()).unwrap_or(0.0)
}

fn read_species(path: &str) -> Option<SpeciesReport> {
    if path.is_empty() || !Path::new(path).exists() {
        return None;
    }
    let content = std::fs::read_to_string(path).ok()?;
    let lines: Vec<&str> = content.lines().filter(|l| !l.trim().is_empty()).collect();
    if lines.is_empty() {
        return None;
    }
    let head: Vec<&str> = lines[0].split('\t').collect();
    let col = |n: &str| head.iter().position(|h| *h == n);
    let get = |f: &[&str], i: Option<usize>| -> String {
        i.and_then(|i| f.get(i)).map(|s| s.to_string()).unwrap_or_default()
    };
    let (i_rank, i_taxid, i_name) = (col("rank"), col("taxid"), col("name"));
    let (i_obs, i_exp, i_enr) = (col("observed"), col("expected"), col("enrichment"));
    let (i_logp, i_q) = (col("log_pvalue"), col("qvalue"));

    let mut taxa: Vec<Taxon> = Vec::new();
    for line in &lines[1..] {
        let f: Vec<&str> = line.split('\t').collect();
        taxa.push(Taxon {
            rank: get(&f, i_rank),
            taxid: get(&f, i_taxid).parse().unwrap_or(0),
            name: get(&f, i_name),
            observed: num(&get(&f, i_obs)),
            expected: num(&get(&f, i_exp)),
            enrichment: num(&get(&f, i_enr)),
            log_p: num(&get(&f, i_logp)),
            q: get(&f, i_q).parse().ok().filter(|v: &f64| v.is_finite()).unwrap_or(1.0),
        });
    }

    // Show only taxa with actual hits: the CLI writes a row for every reference
    // taxon, but a taxon with zero observed k-mers is noise in the report.
    taxa.retain(|t| t.observed > 0.0);

    // Rank by significance (log p-value), not enrichment; ties broken by the
    // count. See the note in species.ts.
    taxa.sort_by(|a, b| {
        a.log_p
            .partial_cmp(&b.log_p)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.observed.partial_cmp(&a.observed).unwrap_or(std::cmp::Ordering::Equal))
    });

    let empty = taxa.is_empty();
    Some(SpeciesReport { path: path.to_string(), taxa, empty })
}

// Layout: "FTXI"/"FTX2" | u32 version | u32 k | ... k-mer count. v1 (FTXI) put
// n at offset 12; v2 (FTX2) put n_taxa at 12 and n_kmers at 16. A real index is
// many MB, so reading 24 bytes always succeeds unless the file is truncated.
fn read_taxdb_info(path: &str) -> Option<TaxdbInfo> {
    if path.is_empty() || !Path::new(path).exists() {
        return None;
    }
    let mut f = std::fs::File::open(path).ok()?;
    let mut buf = [0u8; 24];
    f.read_exact(&mut buf).ok()?;
    let magic = &buf[0..4];
    let k = u32::from_le_bytes(buf[8..12].try_into().ok()?);
    let kmers = if magic == b"FTX2" {
        u64::from_le_bytes(buf[16..24].try_into().ok()?)
    } else if magic == b"FTXI" {
        u64::from_le_bytes(buf[12..20].try_into().ok()?)
    } else {
        return None;
    };
    Some(TaxdbInfo { path: path.to_string(), k, kmers })
}

// The same search order the CLI uses (taxonomyDir_ in FASTag.cpp), so the GUI
// reports on the index the run will actually load.
fn bundled_taxdb(app: &AppHandle) -> Option<String> {
    let bin = crate::fastag::resolve_binary(app).bin;
    let bin_dir = bin.parent()?;
    let mut dirs: Vec<PathBuf> = Vec::new();
    if let Ok(env) = std::env::var("FASTAG_TAXONOMY_DIR") {
        if !env.is_empty() {
            dirs.push(PathBuf::from(env));
        }
    }
    dirs.push(bin_dir.join("..").join("share").join("FASTag").join("taxonomy"));
    dirs.push(bin_dir.join("share-FASTag-taxonomy"));
    dirs.push(bin_dir.join("taxonomy"));
    for d in dirs {
        let p = d.join("tax_k7.taxdb");
        if p.exists() {
            return Some(p.display().to_string());
        }
    }
    None
}

#[tauri::command]
pub fn species(path: String) -> Option<SpeciesReport> {
    read_species(&path)
}

#[tauri::command]
pub fn taxdb_info(app: AppHandle, explicit: Option<String>) -> Option<TaxdbInfo> {
    let path = match explicit {
        Some(e) if !e.trim().is_empty() => e,
        _ => bundled_taxdb(&app)?,
    };
    read_taxdb_info(&path)
}
