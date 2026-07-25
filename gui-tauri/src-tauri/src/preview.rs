// Bounded TSV preview for the results pane. A full run is 100k-1M+ rows; this
// streams the header plus the first N data rows and stops early, so a huge file
// never loads into memory. Ported from gui/src/main/preview.ts.

use serde::Serialize;
use std::io::{BufRead, BufReader};

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Preview {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
    truncated: bool,
    shown: usize,
}

#[tauri::command]
pub fn preview(path: String, max_rows: Option<usize>) -> Preview {
    let max = max_rows.unwrap_or(200);
    let mut header: Vec<String> = Vec::new();
    let mut rows: Vec<Vec<String>> = Vec::new();
    let mut truncated = false;

    if let Ok(file) = std::fs::File::open(&path) {
        let reader = BufReader::new(file);
        let mut first = true;
        for line in reader.lines() {
            let Ok(line) = line else { break };
            if first {
                header = line.split('\t').map(|s| s.to_string()).collect();
                first = false;
                continue;
            }
            if rows.len() < max {
                rows.push(line.split('\t').map(|s| s.to_string()).collect());
            } else {
                truncated = true;
                break;
            }
        }
    }

    let shown = rows.len();
    Preview { header, rows, truncated, shown }
}
