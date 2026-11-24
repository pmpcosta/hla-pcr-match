#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod app;
mod download;
mod export;
mod fasta;
mod iupac;
mod matcher;

use eframe::egui;

fn main() -> eframe::Result {
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1240.0, 780.0])
            .with_min_inner_size([860.0, 560.0])
            .with_title("HLA PCR Match — IPD-IMGT/HLA In-Silico PCR"),
        ..Default::default()
    };

    eframe::run_native(
        "HLA PCR Match",
        native_options,
        Box::new(|cc| Ok(Box::new(app::HlaPrimerApp::new(cc)))),
    )
}
