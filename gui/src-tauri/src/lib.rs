mod fastag;
mod preview;
mod settings;
mod species;

use fastag::RunManager;
use tauri::{Manager, WindowEvent};

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        // single-instance must be registered first: a second launch focuses the
        // existing window instead of racing over the settings file and run state.
        .plugin(tauri_plugin_single_instance::init(|app, _args, _cwd| {
            if let Some(w) = app.get_webview_window("main") {
                let _ = w.unminimize();
                let _ = w.set_focus();
            }
        }))
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .manage(RunManager::default())
        .invoke_handler(tauri::generate_handler![
            fastag::probe,
            fastag::run,
            fastag::cancel,
            preview::preview,
            species::species,
            species::taxdb_info,
            settings::load_settings,
            settings::save_last,
            settings::save_preset,
            settings::delete_preset,
        ])
        .on_window_event(|window, event| {
            // Closing the window must not orphan a running CLI process.
            if let WindowEvent::Destroyed = event {
                let app = window.app_handle().clone();
                let child = app
                    .state::<RunManager>()
                    .current
                    .lock()
                    .unwrap()
                    .as_ref()
                    .map(|r| r.child.clone());
                if let Some(child) = child {
                    if let Ok(mut c) = child.lock() {
                        let _ = c.kill();
                    }
                }
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
