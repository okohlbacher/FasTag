// FasTag CLI integration: locate the binary, run it, stream its stderr as
// events, and cancel by killing the child. Ported from the Electron main
// process (gui/src/main/fastag.ts); the frontend contract is identical.
//
// Security: the binary is spawned with an explicit argv array (never a shell
// string), and build_args() only ever emits flags the tool declares in its
// generated manifest — a crafted invoke cannot inject an unknown option.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use tauri::{AppHandle, Emitter, Manager, State};

// name -> declared type, straight from the tool's own -write_ini. This is the
// allowlist build_args() validates against. Embedded at compile time from the
// single source of truth the frontend also reads.
const MANIFEST: &str = include_str!("../../src/params.generated.json");

fn exe_name() -> &'static str {
    if cfg!(windows) {
        "FasTag.exe"
    } else {
        "FasTag"
    }
}

// The single in-flight run. The batch queue runs sequentially, so at most one
// child exists at a time; a second `run` while one is live is refused.
#[derive(Default)]
pub struct RunManager {
    pub current: Mutex<Option<CurrentRun>>,
    pub seq: AtomicU64,
}

pub struct CurrentRun {
    pub run_id: u64,
    pub child: Arc<Mutex<Child>>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BinaryInfo {
    bin: String,
    data_path: Option<String>,
    source: String,
    ok: bool,
    version: Option<String>,
    detail: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunStarted {
    started: bool,
    run_id: Option<u64>,
    reason: Option<String>,
}

#[derive(Deserialize)]
pub struct RunParams {
    #[serde(rename = "in")]
    input: String,
    out: String,
    #[serde(default)]
    params: Map<String, Value>,
}

pub struct Resolved {
    pub bin: PathBuf,
    pub data: Option<PathBuf>,
    // The bundled taxonomy dir (share/FasTag/taxonomy), if present. Passed to the
    // child as FASTAG_TAXONOMY_DIR so -species works without depending on the
    // CLI's own <bin>/../share path arithmetic (which a copied/symlinked binary
    // can get wrong).
    pub taxonomy: Option<PathBuf>,
    pub source: &'static str,
}

// Resolve the binary in priority order: FASTAG_BIN env override, then bundled
// beside the app (prod resource dir, or the dev resources dir baked in at
// compile time), then the bare name on PATH.
pub fn resolve_binary(app: &AppHandle) -> Resolved {
    if let Ok(env_bin) = std::env::var("FASTAG_BIN") {
        if !env_bin.is_empty() && PathBuf::from(&env_bin).exists() {
            let data = std::env::var("OPENMS_DATA_PATH").ok().map(PathBuf::from);
            // The CLI already honours a user-set FASTAG_TAXONOMY_DIR itself.
            return Resolved { bin: PathBuf::from(env_bin), data, taxonomy: None, source: "env" };
        }
    }
    let mut roots: Vec<PathBuf> = Vec::new();
    if let Ok(rd) = app.path().resource_dir() {
        roots.push(rd.join("resources").join("fastag"));
    }
    // Dev fallback: the resources dir next to this crate (tauri dev does not
    // populate resource_dir the way a bundle does).
    roots.push(PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/resources/fastag")));
    for root in roots {
        let bin = root.join("bin").join(exe_name());
        if bin.exists() {
            let data = root.join("share").join("OpenMS");
            let data = if data.exists() { Some(data) } else { None };
            let taxonomy = root.join("share").join("FasTag").join("taxonomy");
            let taxonomy = if taxonomy.join("tax_k7.taxdb").exists() { Some(taxonomy) } else { None };
            return Resolved { bin, data, taxonomy, source: "bundled" };
        }
    }
    Resolved { bin: PathBuf::from(exe_name()), data: None, taxonomy: None, source: "path" }
}

// TOPP --help prints "... Version: X.Y.Z ...". Capture the version token as
// proof the binary actually runs and links its libraries on this machine.
fn parse_version(text: &str) -> Option<String> {
    let idx = text.find("Version:")?;
    let rest = text[idx + "Version:".len()..].trim_start();
    let tok: String = rest.chars().take_while(|c| !c.is_whitespace() && *c != ',').collect();
    if tok.chars().next().map_or(false, |c| c.is_ascii_digit()) {
        Some(tok)
    } else {
        None
    }
}

#[tauri::command]
pub fn probe(app: AppHandle) -> BinaryInfo {
    let r = resolve_binary(&app);
    let bin_s = r.bin.display().to_string();
    let data_s = r.data.as_ref().map(|d| d.display().to_string());
    let mut cmd = Command::new(&r.bin);
    cmd.arg("--help");
    if let Some(d) = &r.data {
        cmd.env("OPENMS_DATA_PATH", d);
    }
    match cmd.output() {
        Ok(out) => {
            let text = format!(
                "{}\n{}",
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr)
            );
            let version = parse_version(&text);
            let detail = match &version {
                Some(v) => format!("FasTag {v} ({})", r.source),
                None => format!("runs, version unknown ({})", r.source),
            };
            BinaryInfo { bin: bin_s, data_path: data_s, source: r.source.into(), ok: true, version, detail }
        }
        Err(e) => BinaryInfo {
            bin: bin_s,
            data_path: data_s,
            source: r.source.into(),
            ok: false,
            version: None,
            detail: format!("cannot execute: {e}"),
        },
    }
}

// In the manifest but not settable on the command line, or managed by the app.
const NOT_SETTABLE: [&str; 8] = ["version", "log", "debug", "no_progress", "force", "test", "in", "out"];

fn value_to_scalar(v: &Value) -> Option<String> {
    match v {
        Value::String(s) => Some(s.clone()),
        Value::Number(n) => Some(n.to_string()),
        Value::Bool(b) => Some(b.to_string()),
        _ => None, // null / array / object have no scalar CLI form
    }
}

// Build the argv from typed params. Only names the tool declares are emitted,
// arrays only for list-typed options, and each value is its own argv element.
fn build_args(p: &RunParams) -> Vec<String> {
    let manifest: Value = serde_json::from_str(MANIFEST).unwrap_or(Value::Null);
    let mut types: HashMap<String, String> = HashMap::new();
    if let Some(arr) = manifest.get("params").and_then(|v| v.as_array()) {
        for it in arr {
            if let (Some(n), Some(t)) = (
                it.get("name").and_then(|v| v.as_str()),
                it.get("type").and_then(|v| v.as_str()),
            ) {
                types.insert(n.to_string(), t.to_string());
            }
        }
    }
    let not_settable: HashSet<&str> = NOT_SETTABLE.into_iter().collect();

    let mut args: Vec<String> =
        vec!["-in".into(), p.input.clone(), "-out".into(), p.out.clone()];

    for (name, value) in &p.params {
        let Some(t) = types.get(name) else { continue };
        if not_settable.contains(name.as_str()) {
            continue;
        }
        if t == "bool" {
            // presence-only flags: passing "false" would set them
            if value.as_bool() == Some(true) {
                args.push(format!("-{name}"));
            }
            continue;
        }
        if let Some(arr) = value.as_array() {
            // multiple values are only meaningful for a list-typed option
            if t != "string-list" {
                continue;
            }
            let items: Vec<String> = arr
                .iter()
                .filter_map(value_to_scalar)
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();
            if !items.is_empty() {
                args.push(format!("-{name}"));
                args.extend(items);
            }
            continue;
        }
        let Some(s) = value_to_scalar(value) else { continue };
        let s = s.trim().to_string();
        if s.is_empty() {
            continue; // unset optional (an empty path is not "no path")
        }
        args.push(format!("-{name}"));
        args.push(s);
    }

    // The GUI always wants machine-readable progress.
    args.push("-progress".into());
    args
}

const PROGRESS_PREFIX: &str = "FASTAG_PROGRESS ";

fn parse_progress_line(line: &str) -> Option<(u64, u64)> {
    let rest = line.strip_prefix(PROGRESS_PREFIX)?;
    let mut done = None;
    let mut total = None;
    for tok in rest.split_whitespace() {
        if let Some(v) = tok.strip_prefix("done=") {
            done = v.parse().ok();
        } else if let Some(v) = tok.strip_prefix("total=") {
            total = v.parse().ok();
        }
    }
    Some((done?, total?))
}

fn read_stream<R: Read>(stream: R, app: &AppHandle, parse_progress: bool) {
    let reader = BufReader::new(stream);
    for line in reader.lines() {
        let Ok(line) = line else { break };
        if parse_progress {
            if let Some((done, total)) = parse_progress_line(&line) {
                let _ = app.emit("fastag:progress", serde_json::json!({ "done": done, "total": total }));
                continue;
            }
        }
        let _ = app.emit("fastag:log", line);
    }
}

#[tauri::command]
pub fn run(app: AppHandle, state: State<'_, RunManager>, params: RunParams) -> RunStarted {
    {
        let cur = state.current.lock().unwrap();
        if cur.is_some() {
            return RunStarted { started: false, run_id: None, reason: Some("a run is already in progress".into()) };
        }
    }
    let r = resolve_binary(&app);
    let args = build_args(&params);
    let mut cmd = Command::new(&r.bin);
    cmd.args(&args).stdout(Stdio::piped()).stderr(Stdio::piped());
    if let Some(d) = &r.data {
        cmd.env("OPENMS_DATA_PATH", d);
    }
    if let Some(t) = &r.taxonomy {
        // Point -species at the bundled taxonomy explicitly (the CLI's own
        // <bin>/../share lookup fails for a copied/symlinked binary).
        cmd.env("FASTAG_TAXONOMY_DIR", t);
    }
    let mut child = match cmd.spawn() {
        Ok(c) => c,
        // A spawn failure is known synchronously: report it as a non-start so
        // the frontend surfaces "could not start" and settles cleanly.
        Err(e) => return RunStarted { started: false, run_id: None, reason: Some(format!("cannot start: {e}")) },
    };

    let run_id = state.seq.fetch_add(1, Ordering::SeqCst) + 1;
    let stdout = child.stdout.take();
    let stderr = child.stderr.take();
    let child_arc = Arc::new(Mutex::new(child));
    *state.current.lock().unwrap() = Some(CurrentRun { run_id, child: child_arc.clone() });

    // One reader thread per stream (stderr carries the progress contract; stdout
    // is drained as log so a full pipe can never block the child).
    let a_err = app.clone();
    let t_err = std::thread::spawn(move || {
        if let Some(s) = stderr {
            read_stream(s, &a_err, true);
        }
    });
    let a_out = app.clone();
    let t_out = std::thread::spawn(move || {
        if let Some(s) = stdout {
            read_stream(s, &a_out, false);
        }
    });

    // Coordinator: once both streams hit EOF the process has ended; reap it,
    // emit exactly one terminal event, and release the slot. Emitting to a
    // closed window is a harmless no-op (unlike Electron's wc.send).
    let a_done = app.clone();
    std::thread::spawn(move || {
        let _ = t_err.join();
        let _ = t_out.join();
        let status = {
            let mut c = child_arc.lock().unwrap();
            c.wait()
        };
        let (ok, code) = match status {
            Ok(s) => (s.success(), s.code()),
            Err(_) => (false, None),
        };
        {
            let st = a_done.state::<RunManager>();
            let mut cur = st.current.lock().unwrap();
            if cur.as_ref().map_or(false, |c| c.run_id == run_id) {
                *cur = None;
            }
        }
        let _ = a_done.emit(
            "fastag:done",
            serde_json::json!({ "runId": run_id, "ok": ok, "code": code, "message": Value::Null }),
        );
    });

    RunStarted { started: true, run_id: Some(run_id), reason: None }
}

#[tauri::command]
pub fn cancel(state: State<'_, RunManager>) -> Value {
    let cur = state.current.lock().unwrap();
    if let Some(run) = cur.as_ref() {
        if let Ok(mut child) = run.child.lock() {
            let _ = child.kill(); // SIGKILL; the coordinator reaps and emits done
        }
        serde_json::json!({ "cancelled": true })
    } else {
        serde_json::json!({ "cancelled": false })
    }
}
