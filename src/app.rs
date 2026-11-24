//! Native egui/eframe GUI for HLA Primer Match.

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::mpsc::{self, Receiver};
use std::thread;

use eframe::egui::{self, Color32, RichText, TextStyle};
use egui_extras::{Column, TableBuilder};

use crate::download::{
    DownloadProgress, ReleaseInfo, SequenceType, cached_fasta_path, download_imgt_hla,
    read_cached_release_info,
};
use crate::export::{
    export_amplicons_fasta, export_coverage_csv, export_hits_delimited, export_missed_csv,
};
use crate::fasta::AlleleDatabase;
use crate::iupac::{bases_match, clean_and_validate_primer, encode_primer, primer_degeneracy, reverse_complement};
use crate::matcher::{MatchParams, MatchResults, run_match};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ActiveTab {
    Matches,
    Coverage,
    Missed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HitSortCol {
    Allele,
    Start,
    End,
    Length,
    FwdMm,
    RevMm,
    TotalMm,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MissSortCol {
    Allele,
    Locus,
    Reason,
    FwdHits,
    RevHits,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DataSourceMode {
    ImgtHla,
    CustomFasta,
}

enum BgEvent {
    DownloadProgress(DownloadProgress),
    DownloadFinished(Result<ReleaseInfo, String>, bool, bool),
    DbLoaded {
        seq_type: Option<SequenceType>,
        result: Result<Arc<AlleleDatabase>, String>,
    },
    MatchFinished(Result<MatchResults, String>),
}

pub struct HlaPrimerApp {
    // Form fields (matching IPD-MHC Primer Match)
    query: String,
    primer_fwd: String,
    primer_rvs: String,
    seq_type: SequenceType,
    end_identity: usize,
    max_mismatch: usize,
    flip_reverse: bool,
    product_size_min: i64,
    product_size_max: i64,

    // Database state
    source_mode: DataSourceMode,
    release_info: ReleaseInfo,
    cds_db: Option<Arc<AlleleDatabase>>,
    genomic_db: Option<Arc<AlleleDatabase>>,
    custom_db: Option<Arc<AlleleDatabase>>,
    loading_db: bool,

    // Background worker state
    bg_rx: Option<Receiver<BgEvent>>,
    download_progress: Option<DownloadProgress>,
    matching_in_progress: bool,

    // Results & view state
    results: Option<MatchResults>,
    active_tab: ActiveTab,
    filter_text: String,
    selected_hit_idx: Option<usize>,
    hit_sort_col: HitSortCol,
    hit_sort_asc: bool,
    miss_sort_col: MissSortCol,
    miss_sort_asc: bool,
    filtered_hit_indices: Vec<usize>,
    filtered_miss_indices: Vec<usize>,
    indices_dirty: bool,

    // Feedback banners
    error_msg: Option<String>,
    status_msg: Option<String>,
    dark_mode: bool,
}

impl HlaPrimerApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let dark_mode = true;
        cc.egui_ctx.set_visuals(egui::Visuals::dark());

        let release_info = read_cached_release_info();
        let mut app = Self {
            query: String::new(),
            primer_fwd: String::new(),
            primer_rvs: String::new(),
            seq_type: SequenceType::Cds,
            end_identity: 3,
            max_mismatch: 3,
            flip_reverse: true,
            product_size_min: -1,
            product_size_max: -1,

            source_mode: DataSourceMode::ImgtHla,
            release_info,
            cds_db: None,
            genomic_db: None,
            custom_db: None,
            loading_db: false,

            bg_rx: None,
            download_progress: None,
            matching_in_progress: false,

            results: None,
            active_tab: ActiveTab::Matches,
            filter_text: String::new(),
            selected_hit_idx: None,
            hit_sort_col: HitSortCol::Allele,
            hit_sort_asc: true,
            miss_sort_col: MissSortCol::Allele,
            miss_sort_asc: true,
            filtered_hit_indices: Vec::new(),
            filtered_miss_indices: Vec::new(),
            indices_dirty: false,

            error_msg: None,
            status_msg: None,
            dark_mode,
        };

        // Auto-load cached CDS database if present
        let cds_path = cached_fasta_path(SequenceType::Cds);
        if cds_path.exists() {
            app.start_load_cached_db(SequenceType::Cds, cc.egui_ctx.clone());
        } else {
            let gen_path = cached_fasta_path(SequenceType::Genomic);
            if gen_path.exists() {
                app.seq_type = SequenceType::Genomic;
                app.start_load_cached_db(SequenceType::Genomic, cc.egui_ctx.clone());
            }
        }

        app
    }

    fn active_db(&self) -> Option<&Arc<AlleleDatabase>> {
        match self.source_mode {
            DataSourceMode::CustomFasta => self.custom_db.as_ref(),
            DataSourceMode::ImgtHla => match self.seq_type {
                SequenceType::Cds => self.cds_db.as_ref(),
                SequenceType::Genomic => self.genomic_db.as_ref(),
            },
        }
    }

    fn is_busy(&self) -> bool {
        self.loading_db || self.download_progress.is_some() || self.matching_in_progress
    }

    fn start_load_cached_db(&mut self, seq_type: SequenceType, ctx: egui::Context) {
        let path = cached_fasta_path(seq_type);
        if !path.exists() {
            return;
        }
        self.loading_db = true;
        self.error_msg = None;
        self.status_msg = Some(format!("Loading {}...", seq_type.label()));

        let (tx, rx) = mpsc::channel();
        self.bg_rx = Some(rx);
        let label = format!("IMGT/HLA {}", seq_type.label());

        thread::spawn(move || {
            let res = AlleleDatabase::load_from_path(&path, label).map(Arc::new);
            let _ = tx.send(BgEvent::DbLoaded {
                seq_type: Some(seq_type),
                result: res,
            });
            ctx.request_repaint();
        });
    }

    fn start_load_custom_fasta(&mut self, path: PathBuf, ctx: egui::Context) {
        self.loading_db = true;
        self.error_msg = None;
        let fname = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("custom.fasta")
            .to_string();
        self.status_msg = Some(format!("Loading {fname}..."));

        let (tx, rx) = mpsc::channel();
        self.bg_rx = Some(rx);

        thread::spawn(move || {
            let res = AlleleDatabase::load_from_path(&path, fname).map(Arc::new);
            let _ = tx.send(BgEvent::DbLoaded {
                seq_type: None,
                result: res,
            });
            ctx.request_repaint();
        });
    }

    fn start_download(&mut self, cds: bool, genomic: bool, ctx: egui::Context) {
        if self.is_busy() {
            return;
        }
        self.error_msg = None;
        self.status_msg = None;
        self.download_progress = Some(DownloadProgress {
            stage: "Connecting to GitHub (ANHIG/IMGTHLA)...".to_string(),
            bytes_done: 0,
            bytes_total: None,
        });

        let (tx, rx) = mpsc::channel();
        self.bg_rx = Some(rx);

        thread::spawn(move || {
            let (prog_tx, prog_rx) = mpsc::channel::<DownloadProgress>();
            let tx_fwd = tx.clone();
            let ctx_fwd = ctx.clone();

            let forwarder = thread::spawn(move || {
                while let Ok(p) = prog_rx.recv() {
                    let _ = tx_fwd.send(BgEvent::DownloadProgress(p));
                    ctx_fwd.request_repaint();
                }
            });

            let result = download_imgt_hla(cds, genomic, &prog_tx);
            drop(prog_tx);
            let _ = forwarder.join();
            let _ = tx.send(BgEvent::DownloadFinished(result, cds, genomic));
            ctx.request_repaint();
        });
    }

    fn start_match(&mut self, ctx: egui::Context) {
        if self.is_busy() {
            return;
        }
        let Some(db) = self.active_db().cloned() else {
            self.error_msg = Some(
                "No database loaded for the selected Type. Download IMGT/HLA or open a FASTA file."
                    .to_string(),
            );
            return;
        };

        let params = MatchParams {
            query: self.query.clone(),
            primer_fwd: self.primer_fwd.clone(),
            primer_rvs: self.primer_rvs.clone(),
            end_identity: self.end_identity,
            max_mismatch: self.max_mismatch,
            flip_reverse: self.flip_reverse,
            product_size_min: self.product_size_min,
            product_size_max: self.product_size_max,
        };

        self.matching_in_progress = true;
        self.error_msg = None;
        self.status_msg = Some("Running primer match...".to_string());

        let (tx, rx) = mpsc::channel();
        self.bg_rx = Some(rx);

        thread::spawn(move || {
            let res = run_match(&db, &params);
            let _ = tx.send(BgEvent::MatchFinished(res));
            ctx.request_repaint();
        });
    }

    fn poll_background(&mut self, ctx: &egui::Context) {
        let mut events = Vec::new();
        let mut disconnected = false;
        if let Some(ref rx) = self.bg_rx {
            loop {
                match rx.try_recv() {
                    Ok(ev) => events.push(ev),
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => {
                        disconnected = true;
                        break;
                    }
                }
            }
        }
        if disconnected {
            self.bg_rx = None;
        }

        for ev in events {
            match ev {
                BgEvent::DownloadProgress(p) => {
                    self.download_progress = Some(p);
                }
                BgEvent::DownloadFinished(res, downloaded_cds, downloaded_gen) => {
                    self.download_progress = None;
                    match res {
                        Ok(info) => {
                            self.release_info = info;
                            self.source_mode = DataSourceMode::ImgtHla;
                            if downloaded_cds {
                                self.cds_db = None;
                            }
                            if downloaded_gen {
                                self.genomic_db = None;
                            }
                            // Load whichever matches the current dropdown (or the one downloaded)
                            let to_load = match self.seq_type {
                                SequenceType::Cds if downloaded_cds => SequenceType::Cds,
                                SequenceType::Genomic if downloaded_gen => SequenceType::Genomic,
                                _ if downloaded_cds => {
                                    self.seq_type = SequenceType::Cds;
                                    SequenceType::Cds
                                }
                                _ => {
                                    self.seq_type = SequenceType::Genomic;
                                    SequenceType::Genomic
                                }
                            };
                            self.start_load_cached_db(to_load, ctx.clone());
                        }
                        Err(err) => {
                            self.error_msg = Some(format!("Download failed: {err}"));
                        }
                    }
                }
                BgEvent::DbLoaded { seq_type, result } => {
                    self.loading_db = false;
                    match result {
                        Ok(db) => {
                            let count = db.len();
                            let loci = db.loci.len();
                            let label = db.source_label.clone();
                            match seq_type {
                                Some(SequenceType::Cds) => {
                                    self.cds_db = Some(db);
                                    self.source_mode = DataSourceMode::ImgtHla;
                                }
                                Some(SequenceType::Genomic) => {
                                    self.genomic_db = Some(db);
                                    self.source_mode = DataSourceMode::ImgtHla;
                                }
                                None => {
                                    self.custom_db = Some(db);
                                    self.source_mode = DataSourceMode::CustomFasta;
                                }
                            }
                            self.status_msg = Some(format!(
                                "Loaded {label}: {count} alleles across {loci} loci."
                            ));
                        }
                        Err(err) => {
                            self.error_msg = Some(err);
                            self.status_msg = None;
                        }
                    }
                }
                BgEvent::MatchFinished(res) => {
                    self.matching_in_progress = false;
                    match res {
                        Ok(match_res) => {
                            self.status_msg = Some(format!(
                                "Matched {}/{} alleles ({} amplicons) in {} ms.",
                                match_res.matched_unique_alleles,
                                match_res.total_searched,
                                match_res.hits.len(),
                                match_res.elapsed_ms
                            ));
                            self.results = Some(match_res);
                            self.selected_hit_idx = None;
                            self.indices_dirty = true;
                        }
                        Err(err) => {
                            self.error_msg = Some(err);
                            self.status_msg = None;
                        }
                    }
                }
            }
        }
    }

    fn rebuild_sorted_indices(&mut self) {
        let Some(ref res) = self.results else {
            self.filtered_hit_indices.clear();
            self.filtered_miss_indices.clear();
            self.indices_dirty = false;
            return;
        };

        let needle = self.filter_text.trim().to_ascii_lowercase();

        self.filtered_hit_indices = res
            .hits
            .iter()
            .enumerate()
            .filter(|(_, h)| {
                needle.is_empty()
                    || h.allele.to_ascii_lowercase().contains(&needle)
                    || h.locus.to_ascii_lowercase().contains(&needle)
                    || h.accession.to_ascii_lowercase().contains(&needle)
            })
            .map(|(i, _)| i)
            .collect();

        let col = self.hit_sort_col;
        let asc = self.hit_sort_asc;
        self.filtered_hit_indices.sort_by(|&a, &b| {
            let ha = &res.hits[a];
            let hb = &res.hits[b];
            let ord = match col {
                HitSortCol::Allele => ha.allele.cmp(&hb.allele),
                HitSortCol::Start => ha.start.cmp(&hb.start),
                HitSortCol::End => ha.end.cmp(&hb.end),
                HitSortCol::Length => ha.length.cmp(&hb.length),
                HitSortCol::FwdMm => ha.fwd_mismatches.cmp(&hb.fwd_mismatches),
                HitSortCol::RevMm => ha.rev_mismatches.cmp(&hb.rev_mismatches),
                HitSortCol::TotalMm => (ha.fwd_mismatches + ha.rev_mismatches)
                    .cmp(&(hb.fwd_mismatches + hb.rev_mismatches)),
            };
            if asc { ord } else { ord.reverse() }
        });

        self.filtered_miss_indices = res
            .missed
            .iter()
            .enumerate()
            .filter(|(_, m)| {
                needle.is_empty()
                    || m.allele.to_ascii_lowercase().contains(&needle)
                    || m.locus.to_ascii_lowercase().contains(&needle)
                    || m.reason.label().to_ascii_lowercase().contains(&needle)
            })
            .map(|(i, _)| i)
            .collect();

        let mcol = self.miss_sort_col;
        let masc = self.miss_sort_asc;
        self.filtered_miss_indices.sort_by(|&a, &b| {
            let ma = &res.missed[a];
            let mb = &res.missed[b];
            let ord = match mcol {
                MissSortCol::Allele => ma.allele.cmp(&mb.allele),
                MissSortCol::Locus => ma.locus.cmp(&mb.locus),
                MissSortCol::Reason => ma.reason.label().cmp(&mb.reason.label()),
                MissSortCol::FwdHits => ma.fwd_hits.cmp(&mb.fwd_hits),
                MissSortCol::RevHits => ma.rev_hits.cmp(&mb.rev_hits),
            };
            if masc { ord } else { ord.reverse() }
        });

        self.indices_dirty = false;
    }

    fn load_example(&mut self) {
        self.query = "HLA-A".to_string();
        // Real conserved HLA-A exon 1 -> exon 2 CDS primers from A*01:01:01:01
        self.primer_fwd = "ATGGCCGTCATGGCGCCCCG".to_string();
        self.primer_rvs = "CTGGTCCCAATACTCCGG".to_string();
        self.seq_type = SequenceType::Cds;
        self.end_identity = 3;
        self.max_mismatch = 3;
        self.flip_reverse = true;
        self.product_size_min = -1;
        self.product_size_max = -1;
        self.error_msg = None;
    }

    fn reset_form(&mut self) {
        self.query.clear();
        self.primer_fwd.clear();
        self.primer_rvs.clear();
        self.end_identity = 3;
        self.max_mismatch = 3;
        self.flip_reverse = true;
        self.product_size_min = -1;
        self.product_size_max = -1;
        self.error_msg = None;
    }

    fn render_top_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("top_bar").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("HLA PCR Match").strong().size(18.0));
                ui.separator();

                // Database status badge
                if let Some(db) = self.active_db() {
                    let badge = match self.source_mode {
                        DataSourceMode::ImgtHla => format!(
                            "{} • {} alleles • {}",
                            self.seq_type.label(),
                            format_num(db.len()),
                            self.release_info.display()
                        ),
                        DataSourceMode::CustomFasta => {
                            format!(
                                "Custom: {} • {} alleles",
                                db.source_label,
                                format_num(db.len())
                            )
                        }
                    };
                    ui.colored_label(Color32::from_rgb(76, 195, 138), badge);
                } else {
                    ui.colored_label(
                        Color32::from_rgb(245, 179, 66),
                        format!("{} not loaded", self.seq_type.label()),
                    );
                }

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let theme_label = if self.dark_mode { "Light" } else { "Dark" };
                    if ui.button(theme_label).clicked() {
                        self.dark_mode = !self.dark_mode;
                        ctx.set_visuals(if self.dark_mode {
                            egui::Visuals::dark()
                        } else {
                            egui::Visuals::light()
                        });
                    }

                    ui.add_enabled_ui(!self.is_busy(), |ui| {
                        if ui.button("Open FASTA...").clicked() {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("FASTA / ZIP", &["fasta", "fa", "fna", "zip"])
                                .pick_file()
                            {
                                self.start_load_custom_fasta(path, ctx.clone());
                            }
                        }

                        ui.menu_button("Download IMGT/HLA ▾", |ui| {
                            if ui.button("CDS + Genomic (~51 MB)").clicked() {
                                ui.close_menu();
                                self.start_download(true, true, ctx.clone());
                            }
                            if ui.button("CDS only — hla_nuc.fasta (~43 MB)").clicked() {
                                ui.close_menu();
                                self.start_download(true, false, ctx.clone());
                            }
                            if ui
                                .button("Genomic only — hla_gen.fasta.zip (~8.4 MB)")
                                .clicked()
                            {
                                ui.close_menu();
                                self.start_download(false, true, ctx.clone());
                            }
                        });
                    });
                });
            });

            // Progress bar during download or loading
            if let Some(ref prog) = self.download_progress {
                ui.add_space(4.0);
                let frac = prog
                    .bytes_total
                    .map(|t| if t > 0 { prog.bytes_done as f32 / t as f32 } else { 0.0 });
                let text = match prog.bytes_total {
                    Some(total) => format!(
                        "{} ({:.1} / {:.1} MB)",
                        prog.stage,
                        prog.bytes_done as f64 / 1_048_576.0,
                        total as f64 / 1_048_576.0
                    ),
                    None => prog.stage.clone(),
                };
                if let Some(f) = frac {
                    ui.add(egui::ProgressBar::new(f).text(text).animate(true));
                } else {
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label(text);
                    });
                }
            } else if self.loading_db || self.matching_in_progress {
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(self.status_msg.as_deref().unwrap_or("Working..."));
                });
            }
            ui.add_space(4.0);
        });
    }

    fn render_form_panel(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("form_panel")
            .default_width(370.0)
            .min_width(310.0)
            .max_width(520.0)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.add_space(6.0);
                    ui.heading("Primer Match Query");
                    ui.add_space(6.0);

                    // 1. Specific sequences required
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Specific sequences required").strong());
                        help_icon(
                            ui,
                            "Choose specific loci or alleles separated by a comma or newline.\n\
                             Examples: HLA-A, DRB1, B*07:02, A*01, MICA.\n\
                             Leave empty to search all alleles in the loaded database.",
                        );
                    });
                    ui.label(
                        RichText::new("(separated by a new line or a comma; empty = all loci)")
                            .small()
                            .weak(),
                    );
                    ui.add(
                        egui::TextEdit::multiline(&mut self.query)
                            .desired_rows(4)
                            .desired_width(f32::INFINITY)
                            .font(TextStyle::Monospace)
                            .hint_text("HLA-A, HLA-B\nDRB1*15:01"),
                    );

                    // Quick locus buttons
                    ui.horizontal_wrapped(|ui| {
                        ui.label(RichText::new("Quick:").small().weak());
                        for loc in ["HLA-A", "HLA-B", "HLA-C", "DRB1", "DQB1", "DPB1"] {
                            if ui.small_button(loc).clicked() {
                                if self.query.trim().is_empty() {
                                    self.query = loc.to_string();
                                } else {
                                    self.query.push_str(", ");
                                    self.query.push_str(loc);
                                }
                            }
                        }
                        if ui.small_button("All").clicked() {
                            self.query.clear();
                        }
                    });

                    ui.add_space(10.0);

                    // 2. Forward primer
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Forward primer (5'→3'):").strong());
                        help_icon(
                            ui,
                            "5'→3' sequence of the forward primer. IUPAC degenerate bases (R, Y, S, W, K, M, B, D, H, V, N) are supported.",
                        );
                    });
                    ui.add(
                        egui::TextEdit::singleline(&mut self.primer_fwd)
                            .desired_width(f32::INFINITY)
                            .font(TextStyle::Monospace)
                            .hint_text("5'-ATGGCCGTCATGGCGCCCCG-3'"),
                    );
                    render_primer_status(ui, &self.primer_fwd, false);

                    ui.add_space(8.0);

                    // 3. Reverse primer
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Reverse primer:").strong());
                        help_icon(
                            ui,
                            "Sequence of the reverse primer (IUPAC degenerate bases supported).\n\
                             When 'Flip reverse' is Yes (default), enter it in standard 5'→3' oligo order and it will be reverse-complemented automatically.",
                        );
                    });
                    ui.add(
                        egui::TextEdit::singleline(&mut self.primer_rvs)
                            .desired_width(f32::INFINITY)
                            .font(TextStyle::Monospace)
                            .hint_text("5'-CTGGTCCCAATACTCCGG-3'"),
                    );
                    render_primer_status(ui, &self.primer_rvs, self.flip_reverse);

                    ui.add_space(10.0);

                    // 4. Type (CDS database / Genomic database)
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Type:").strong());
                        help_icon(
                            ui,
                            "Select coding DNA (CDS: hla_nuc.fasta) or genomic DNA (Genomic: hla_gen.fasta including introns/UTRs).",
                        );
                        let prev_type = self.seq_type;
                        egui::ComboBox::from_id_salt("seq_type_combo")
                            .selected_text(self.seq_type.label())
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut self.seq_type,
                                    SequenceType::Cds,
                                    SequenceType::Cds.label(),
                                );
                                ui.selectable_value(
                                    &mut self.seq_type,
                                    SequenceType::Genomic,
                                    SequenceType::Genomic.label(),
                                );
                            });

                        if self.seq_type != prev_type {
                            self.source_mode = DataSourceMode::ImgtHla;
                            let already_loaded = match self.seq_type {
                                SequenceType::Cds => self.cds_db.is_some(),
                                SequenceType::Genomic => self.genomic_db.is_some(),
                            };
                            if !already_loaded && cached_fasta_path(self.seq_type).exists() {
                                self.start_load_cached_db(self.seq_type, ctx.clone());
                            }
                        }
                    });

                    if self.source_mode == DataSourceMode::CustomFasta {
                        ui.horizontal(|ui| {
                            ui.colored_label(
                                Color32::from_rgb(130, 180, 255),
                                "Using custom FASTA file.",
                            );
                            if ui.small_button("Switch to IMGT/HLA").clicked() {
                                self.source_mode = DataSourceMode::ImgtHla;
                                if self.active_db().is_none()
                                    && cached_fasta_path(self.seq_type).exists()
                                {
                                    self.start_load_cached_db(self.seq_type, ctx.clone());
                                }
                            }
                        });
                    } else if self.active_db().is_none() && !self.is_busy() {
                        ui.horizontal(|ui| {
                            ui.colored_label(
                                Color32::from_rgb(245, 179, 66),
                                format!("{} not cached locally.", self.seq_type.label()),
                            );
                            if ui.small_button("Download now").clicked() {
                                let download_cds = self.seq_type == SequenceType::Cds;
                                let download_genomic = self.seq_type == SequenceType::Genomic;
                                self.start_download(download_cds, download_genomic, ctx.clone());
                            }
                        });
                    }

                    ui.add_space(10.0);

                    // 5. Options group
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        ui.label(RichText::new("Options").strong());
                        ui.add_space(4.0);
                        egui::Grid::new("options_grid")
                            .num_columns(2)
                            .spacing([12.0, 8.0])
                            .show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.label("End Identity:");
                                    help_icon(
                                        ui,
                                        "Number of bases at the 3' end of each primer that must match identically (0 mismatches allowed).",
                                    );
                                });
                                ui.add(
                                    egui::DragValue::new(&mut self.end_identity)
                                        .range(0..=50)
                                        .suffix(" bp"),
                                );
                                ui.end_row();

                                ui.horizontal(|ui| {
                                    ui.label("Max mismatches:");
                                    help_icon(
                                        ui,
                                        "Maximum total mismatches allowed across the entire primer sequence.",
                                    );
                                });
                                ui.add(egui::DragValue::new(&mut self.max_mismatch).range(0..=20));
                                ui.end_row();

                                ui.horizontal(|ui| {
                                    ui.label("Flip reverse:");
                                    help_icon(
                                        ui,
                                        "If the reverse primer is provided in 5'→3' direction, Yes automatically reverse-complements it.",
                                    );
                                });
                                egui::ComboBox::from_id_salt("flip_reverse_combo")
                                    .selected_text(if self.flip_reverse { "Yes" } else { "No" })
                                    .show_ui(ui, |ui| {
                                        ui.selectable_value(&mut self.flip_reverse, true, "Yes");
                                        ui.selectable_value(&mut self.flip_reverse, false, "No");
                                    });
                                ui.end_row();
                            });
                    });

                    ui.add_space(8.0);

                    // 6. Product Size group
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label(RichText::new("Product Size").strong());
                            help_icon(
                                ui,
                                "Predicted amplicon size range in bp. Use -1 for no minimum or maximum limit.",
                            );
                        });
                        ui.add_space(4.0);
                        ui.horizontal(|ui| {
                            ui.label("Min:");
                            ui.add(
                                egui::DragValue::new(&mut self.product_size_min)
                                    .range(-1..=100_000)
                                    .speed(10),
                            );
                            ui.add_space(12.0);
                            ui.label("Max:");
                            ui.add(
                                egui::DragValue::new(&mut self.product_size_max)
                                    .range(-1..=100_000)
                                    .speed(10),
                            );
                        });
                        ui.label(RichText::new("(-1 = unbounded)").small().weak());
                    });

                    ui.add_space(12.0);

                    // Action buttons
                    let ctrl_enter =
                        ctx.input(|i| i.modifiers.ctrl && i.key_pressed(egui::Key::Enter));
                    ui.horizontal(|ui| {
                        let submit_btn = egui::Button::new(
                            RichText::new("  Submit  ").strong().color(Color32::WHITE),
                        )
                        .fill(Color32::from_rgb(38, 128, 214))
                        .min_size(egui::vec2(110.0, 32.0));

                        let clicked = ui.add_enabled(!self.is_busy(), submit_btn).clicked();
                        if clicked || (ctrl_enter && !self.is_busy()) {
                            self.start_match(ctx.clone());
                        }

                        if ui
                            .add_sized([76.0, 32.0], egui::Button::new("Example"))
                            .clicked()
                        {
                            self.load_example();
                        }

                        if ui
                            .add_sized([68.0, 32.0], egui::Button::new("Reset"))
                            .clicked()
                        {
                            self.reset_form();
                        }
                    });

                    if let Some(ref err) = self.error_msg {
                        ui.add_space(10.0);
                        egui::Frame::group(ui.style())
                            .fill(Color32::from_rgba_unmultiplied(200, 50, 50, 35))
                            .show(ui, |ui| {
                                ui.colored_label(Color32::from_rgb(255, 115, 115), err);
                            });
                    }
                });
            });
    }

    fn render_results_panel(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            // Welcome onboarding if no database is loaded and no results exist
            if self.active_db().is_none() && self.results.is_none() {
                ui.vertical_centered(|ui| {
                    ui.add_space(60.0);
                    ui.heading("Welcome to HLA PCR Match");
                    ui.add_space(8.0);
                    ui.label(
                        "In-silico PCR primer matching against the official IPD-IMGT/HLA human allele database.",
                    );
                    ui.add_space(16.0);
                    ui.add_enabled_ui(!self.is_busy(), |ui| {
                        let btn = egui::Button::new(
                            RichText::new("Download Latest IPD-IMGT/HLA (CDS + Genomic)")
                                .strong()
                                .color(Color32::WHITE),
                        )
                        .fill(Color32::from_rgb(38, 128, 214))
                        .min_size(egui::vec2(320.0, 40.0));
                        if ui.add(btn).clicked() {
                            self.start_download(true, true, ctx.clone());
                        }
                        ui.add_space(8.0);
                        if ui.button("Download CDS Only (faster, ~43 MB)").clicked() {
                            self.start_download(true, false, ctx.clone());
                        }
                        ui.add_space(4.0);
                        if ui.button("Or Open a Local FASTA / ZIP File...").clicked() {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("FASTA / ZIP", &["fasta", "fa", "fna", "zip"])
                                .pick_file()
                            {
                                self.start_load_custom_fasta(path, ctx.clone());
                            }
                        }
                    });
                });
                return;
            }

            if self.results.is_none() {
                ui.vertical_centered(|ui| {
                    ui.add_space(60.0);
                    ui.heading("Ready to Match");
                    ui.add_space(8.0);
                    ui.label(
                        "Enter a forward and reverse primer on the left (or click Example), then press Submit.",
                    );
                    if let Some(db) = self.active_db() {
                        ui.add_space(12.0);
                        ui.label(
                            RichText::new(format!(
                                "Loaded {} alleles across {} loci: {}",
                                format_num(db.len()),
                                db.loci.len(),
                                db.loci.join(", ")
                            ))
                            .small()
                            .weak(),
                        );
                    }
                });
                return;
            }

            if self.indices_dirty {
                self.rebuild_sorted_indices();
            }

            // Summary metrics & export bar
            let (matched_u, total_s, hit_cnt, miss_cnt, locus_cnt, elapsed, unmatched_q) = {
                let res = self.results.as_ref().unwrap();
                (
                    res.matched_unique_alleles,
                    res.total_searched,
                    res.hits.len(),
                    res.missed.len(),
                    res.locus_coverage.len(),
                    res.elapsed_ms,
                    res.unmatched_queries.clone(),
                )
            };
            let cov_pct = if total_s > 0 {
                (matched_u as f64 / total_s as f64) * 100.0
            } else {
                0.0
            };

            egui::Frame::group(ui.style()).show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label(
                        RichText::new(format!(
                            "Coverage: {}/{} alleles ({:.1}%)",
                            format_num(matched_u),
                            format_num(total_s),
                            cov_pct
                        ))
                        .strong()
                        .color(if cov_pct >= 90.0 {
                            Color32::from_rgb(76, 195, 138)
                        } else if cov_pct > 0.0 {
                            Color32::from_rgb(245, 179, 66)
                        } else {
                            Color32::from_rgb(255, 115, 115)
                        }),
                    );
                    ui.separator();
                    ui.label(format!("Amplicons: {}", format_num(hit_cnt)));
                    ui.separator();
                    ui.label(format!("Missed: {}", format_num(miss_cnt)));
                    ui.separator();
                    ui.label(RichText::new(format!("{elapsed} ms")).weak());

                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.menu_button("Export ▾", |ui| {
                            if ui.button("Matched Hits (.csv)").clicked() {
                                ui.close_menu();
                                self.do_export_hits(',');
                            }
                            if ui.button("Matched Hits (.tsv)").clicked() {
                                ui.close_menu();
                                self.do_export_hits('\t');
                            }
                            if ui.button("Missed Alleles (.csv)").clicked() {
                                ui.close_menu();
                                self.do_export_missed();
                            }
                            if ui.button("Locus Coverage (.csv)").clicked() {
                                ui.close_menu();
                                self.do_export_coverage();
                            }
                            if ui.button("Amplicon Sequences (.fasta)").clicked() {
                                ui.close_menu();
                                self.do_export_fasta();
                            }
                        });
                    });
                });

                if !unmatched_q.is_empty() {
                    ui.add_space(4.0);
                    ui.colored_label(
                        Color32::from_rgb(245, 179, 66),
                        format!(
                            "Warning — unrecognized query term(s): {}",
                            unmatched_q.join(", ")
                        ),
                    );
                }
            });

            ui.add_space(6.0);

            // Tabs + search filter
            ui.horizontal(|ui| {
                ui.selectable_value(
                    &mut self.active_tab,
                    ActiveTab::Matches,
                    format!("Matches ({})", format_num(hit_cnt)),
                );
                ui.selectable_value(
                    &mut self.active_tab,
                    ActiveTab::Coverage,
                    format!("Locus Coverage ({locus_cnt})"),
                );
                ui.selectable_value(
                    &mut self.active_tab,
                    ActiveTab::Missed,
                    format!("Missed Alleles ({})", format_num(miss_cnt)),
                );

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if !self.filter_text.is_empty() && ui.small_button("✕").clicked() {
                        self.filter_text.clear();
                        self.indices_dirty = true;
                    }
                    let resp = ui.add(
                        egui::TextEdit::singleline(&mut self.filter_text)
                            .desired_width(210.0)
                            .hint_text("Filter allele / locus..."),
                    );
                    if resp.changed() {
                        self.indices_dirty = true;
                    }
                });
            });

            ui.separator();

            // Reserve bottom inspector space when a hit is selected on the Matches tab
            let show_inspector =
                self.active_tab == ActiveTab::Matches && self.selected_hit_idx.is_some();
            let table_height = if show_inspector {
                (ui.available_height() - 170.0).max(160.0)
            } else {
                ui.available_height().max(160.0)
            };

            ui.allocate_ui_with_layout(
                egui::vec2(ui.available_width(), table_height),
                egui::Layout::top_down(egui::Align::Min),
                |ui| match self.active_tab {
                    ActiveTab::Matches => self.render_matches_table(ui),
                    ActiveTab::Coverage => self.render_coverage_table(ui),
                    ActiveTab::Missed => self.render_missed_table(ui),
                },
            );

            if show_inspector {
                ui.separator();
                self.render_alignment_inspector(ui);
            }
        });
    }

    fn render_matches_table(&mut self, ui: &mut egui::Ui) {
        let Some(ref res) = self.results else {
            return;
        };

        if self.filtered_hit_indices.is_empty() {
            ui.add_space(20.0);
            ui.label("No matching amplicons.");
            return;
        }

        let mut sort_change: Option<HitSortCol> = None;
        let mut newly_selected: Option<usize> = None;

        TableBuilder::new(ui)
            .striped(true)
            .resizable(true)
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
            .sense(egui::Sense::click())
            .column(Column::initial(155.0).at_least(100.0)) // Allele
            .column(Column::initial(65.0).at_least(50.0)) // Start
            .column(Column::initial(65.0).at_least(50.0)) // End
            .column(Column::initial(70.0).at_least(55.0)) // Length
            .column(Column::initial(65.0).at_least(50.0)) // Fwd MM
            .column(Column::initial(65.0).at_least(50.0)) // Rev MM
            .column(Column::initial(70.0).at_least(55.0)) // Total MM
            .column(Column::remainder().at_least(140.0)) // Binding sites preview
            .header(22.0, |mut header| {
                let cols = [
                    ("Allele", HitSortCol::Allele),
                    ("Start", HitSortCol::Start),
                    ("End", HitSortCol::End),
                    ("Length", HitSortCol::Length),
                    ("Fwd MM", HitSortCol::FwdMm),
                    ("Rev MM", HitSortCol::RevMm),
                    ("Total MM", HitSortCol::TotalMm),
                ];
                for (title, col) in cols {
                    header.col(|ui| {
                        let arrow = if self.hit_sort_col == col {
                            if self.hit_sort_asc { " ▲" } else { " ▼" }
                        } else {
                            ""
                        };
                        if ui
                            .button(RichText::new(format!("{title}{arrow}")).strong())
                            .clicked()
                        {
                            sort_change = Some(col);
                        }
                    });
                }
                header.col(|ui| {
                    ui.label(RichText::new("Fwd / Rev Template Sites").strong());
                });
            })
            .body(|body| {
                let row_count = self.filtered_hit_indices.len();
                body.rows(20.0, row_count, |mut row| {
                    let row_i = row.index();
                    let hit_idx = self.filtered_hit_indices[row_i];
                    let hit = &res.hits[hit_idx];
                    let is_selected = self.selected_hit_idx == Some(hit_idx);
                    row.set_selected(is_selected);

                    row.col(|ui| {
                        ui.monospace(&hit.allele);
                    });
                    row.col(|ui| {
                        ui.label(hit.start.to_string());
                    });
                    row.col(|ui| {
                        ui.label(hit.end.to_string());
                    });
                    row.col(|ui| {
                        ui.label(hit.length.to_string());
                    });
                    row.col(|ui| {
                        mm_badge(ui, hit.fwd_mismatches);
                    });
                    row.col(|ui| {
                        mm_badge(ui, hit.rev_mismatches);
                    });
                    row.col(|ui| {
                        mm_badge(ui, hit.fwd_mismatches + hit.rev_mismatches);
                    });
                    row.col(|ui| {
                        ui.monospace(format!("{} … {}", hit.fwd_template, hit.rev_template));
                    });

                    if row.response().clicked() {
                        newly_selected = Some(hit_idx);
                    }
                });
            });

        if let Some(col) = sort_change {
            if self.hit_sort_col == col {
                self.hit_sort_asc = !self.hit_sort_asc;
            } else {
                self.hit_sort_col = col;
                self.hit_sort_asc = true;
            }
            self.indices_dirty = true;
        }
        if let Some(idx) = newly_selected {
            self.selected_hit_idx = Some(idx);
        }
    }

    fn render_coverage_table(&self, ui: &mut egui::Ui) {
        let Some(ref res) = self.results else {
            return;
        };
        let needle = self.filter_text.trim().to_ascii_lowercase();
        let rows: Vec<_> = res
            .locus_coverage
            .iter()
            .filter(|c| needle.is_empty() || c.locus.to_ascii_lowercase().contains(&needle))
            .collect();

        TableBuilder::new(ui)
            .striped(true)
            .resizable(true)
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
            .column(Column::initial(110.0))
            .column(Column::initial(85.0))
            .column(Column::initial(85.0))
            .column(Column::initial(85.0))
            .column(Column::remainder().at_least(180.0))
            .header(22.0, |mut header| {
                for title in ["Locus", "Matched", "Missed", "Total", "Coverage"] {
                    header.col(|ui| {
                        ui.label(RichText::new(title).strong());
                    });
                }
            })
            .body(|body| {
                body.rows(22.0, rows.len(), |mut row| {
                    let c = rows[row.index()];
                    row.col(|ui| {
                        ui.monospace(&c.locus);
                    });
                    row.col(|ui| {
                        ui.label(format_num(c.matched_alleles));
                    });
                    row.col(|ui| {
                        ui.label(format_num(c.missed_alleles));
                    });
                    row.col(|ui| {
                        ui.label(format_num(c.total_alleles));
                    });
                    row.col(|ui| {
                        ui.add(
                            egui::ProgressBar::new((c.coverage_pct / 100.0) as f32)
                                .text(format!("{:.1}%", c.coverage_pct)),
                        );
                    });
                });
            });
    }

    fn render_missed_table(&mut self, ui: &mut egui::Ui) {
        let Some(ref res) = self.results else {
            return;
        };
        if self.filtered_miss_indices.is_empty() {
            ui.add_space(20.0);
            ui.label("No missed alleles (100% coverage in current view).");
            return;
        }

        let mut sort_change: Option<MissSortCol> = None;

        TableBuilder::new(ui)
            .striped(true)
            .resizable(true)
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
            .column(Column::initial(160.0))
            .column(Column::initial(80.0))
            .column(Column::remainder().at_least(220.0))
            .column(Column::initial(80.0))
            .column(Column::initial(80.0))
            .header(22.0, |mut header| {
                let cols = [
                    ("Allele", MissSortCol::Allele),
                    ("Locus", MissSortCol::Locus),
                    ("Diagnostic Reason", MissSortCol::Reason),
                    ("Fwd Sites", MissSortCol::FwdHits),
                    ("Rev Sites", MissSortCol::RevHits),
                ];
                for (title, col) in cols {
                    header.col(|ui| {
                        let arrow = if self.miss_sort_col == col {
                            if self.miss_sort_asc { " ▲" } else { " ▼" }
                        } else {
                            ""
                        };
                        if ui
                            .button(RichText::new(format!("{title}{arrow}")).strong())
                            .clicked()
                        {
                            sort_change = Some(col);
                        }
                    });
                }
            })
            .body(|body| {
                let count = self.filtered_miss_indices.len();
                body.rows(20.0, count, |mut row| {
                    let m = &res.missed[self.filtered_miss_indices[row.index()]];
                    row.col(|ui| {
                        ui.monospace(&m.allele);
                    });
                    row.col(|ui| {
                        ui.monospace(&m.locus);
                    });
                    row.col(|ui| {
                        ui.label(m.reason.label());
                    });
                    row.col(|ui| {
                        ui.label(m.fwd_hits.to_string());
                    });
                    row.col(|ui| {
                        ui.label(m.rev_hits.to_string());
                    });
                });
            });

        if let Some(col) = sort_change {
            if self.miss_sort_col == col {
                self.miss_sort_asc = !self.miss_sort_asc;
            } else {
                self.miss_sort_col = col;
                self.miss_sort_asc = true;
            }
            self.indices_dirty = true;
        }
    }

    fn render_alignment_inspector(&mut self, ui: &mut egui::Ui) {
        let Some(hit_idx) = self.selected_hit_idx else {
            return;
        };
        let Some(ref res) = self.results else {
            return;
        };
        let Some(hit) = res.hits.get(hit_idx) else {
            return;
        };

        let mut close_inspector = false;
        let mut copy_seq: Option<String> = None;

        ui.horizontal(|ui| {
            ui.label(
                RichText::new(format!(
                    "Binding Alignment — {}  ({}..{}, {} bp)",
                    hit.allele, hit.start, hit.end, hit.length
                ))
                .strong(),
            );
            if !hit.accession.is_empty() {
                ui.label(RichText::new(&hit.accession).weak());
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.small_button("Close").clicked() {
                    close_inspector = true;
                }
                if ui.small_button("Copy Amplicon DNA").clicked() {
                    if let Some(db) = self.active_db() {
                        if let Some(rec) = db.records.get(hit.allele_index) {
                            if hit.end <= rec.sequence.len() && hit.start < hit.end {
                                copy_seq = Some(
                                    String::from_utf8_lossy(&rec.sequence[hit.start..hit.end])
                                        .into_owned(),
                                );
                            }
                        }
                    }
                }
            });
        });

        ui.add_space(4.0);
        egui::ScrollArea::horizontal().show(ui, |ui| {
            ui.horizontal(|ui| {
                // Forward alignment block
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    ui.label(
                        RichText::new(format!(
                            "Forward Primer (5'→3') • {} mismatch(es) • green = 3' end identity",
                            hit.fwd_mismatches
                        ))
                        .small(),
                    );
                    render_colored_alignment(
                        ui,
                        &res.fwd_primer_used,
                        &hit.fwd_template,
                        res.end_identity_used,
                        true,
                    );
                });

                ui.add_space(8.0);

                // Reverse alignment block (sense strand orientation)
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    ui.label(
                        RichText::new(format!(
                            "Reverse Motif (sense +) • {} mismatch(es) • green = rev 3' end identity",
                            hit.rev_mismatches
                        ))
                        .small(),
                    );
                    render_colored_alignment(
                        ui,
                        &res.rev_motif_used,
                        &hit.rev_template,
                        res.end_identity_used,
                        false,
                    );
                });
            });
        });

        if let Some(seq) = copy_seq {
            ui.ctx().copy_text(seq);
            self.status_msg = Some("Amplicon sequence copied to clipboard.".to_string());
        }
        if close_inspector {
            self.selected_hit_idx = None;
        }
    }

    fn do_export_hits(&mut self, sep: char) {
        let Some(ref res) = self.results else { return };
        let (ext, name) = if sep == '\t' {
            ("tsv", "hla_pcr_hits.tsv")
        } else {
            ("csv", "hla_pcr_hits.csv")
        };
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name(name)
            .add_filter(ext.to_uppercase(), &[ext])
            .save_file()
        {
            match export_hits_delimited(&path, res, sep) {
                Ok(n) => {
                    self.status_msg = Some(format!("Exported {n} hits to {}", path.display()));
                }
                Err(e) => self.error_msg = Some(e),
            }
        }
    }

    fn do_export_missed(&mut self) {
        let Some(ref res) = self.results else { return };
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("hla_pcr_missed.csv")
            .add_filter("CSV", &["csv"])
            .save_file()
        {
            match export_missed_csv(&path, res) {
                Ok(n) => {
                    self.status_msg =
                        Some(format!("Exported {n} missed alleles to {}", path.display()));
                }
                Err(e) => self.error_msg = Some(e),
            }
        }
    }

    fn do_export_coverage(&mut self) {
        let Some(ref res) = self.results else { return };
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("hla_pcr_coverage.csv")
            .add_filter("CSV", &["csv"])
            .save_file()
        {
            match export_coverage_csv(&path, res) {
                Ok(n) => {
                    self.status_msg =
                        Some(format!("Exported {n} locus rows to {}", path.display()));
                }
                Err(e) => self.error_msg = Some(e),
            }
        }
    }

    fn do_export_fasta(&mut self) {
        let (Some(res), Some(db)) = (&self.results, self.active_db()) else {
            return;
        };
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("hla_pcr_amplicons.fasta")
            .add_filter("FASTA", &["fasta", "fa"])
            .save_file()
        {
            match export_amplicons_fasta(&path, db, res) {
                Ok(n) => {
                    self.status_msg = Some(format!(
                        "Exported {n} amplicon sequences to {}",
                        path.display()
                    ));
                }
                Err(e) => self.error_msg = Some(e),
            }
        }
    }
}

impl eframe::App for HlaPrimerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_background(ctx);
        self.render_top_bar(ctx);
        self.render_form_panel(ctx);
        self.render_results_panel(ctx);

        if let Some(ref msg) = self.status_msg {
            egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(msg).small().weak());
                });
            });
        }
    }
}

fn help_icon(ui: &mut egui::Ui, tooltip: &str) {
    ui.label(RichText::new("(?)").small().weak())
        .on_hover_text(tooltip);
}

fn render_primer_status(ui: &mut egui::Ui, raw: &str, show_rc: bool) {
    if raw.trim().is_empty() {
        return;
    }
    match clean_and_validate_primer(raw, "Primer") {
        Ok(clean) => {
            let deg = primer_degeneracy(&clean);
            let mut info = format!("{} nt • degeneracy {}x", clean.len(), deg);
            if show_rc {
                let rc = reverse_complement(&clean);
                info.push_str(&format!(" • RC: {rc}"));
            }
            ui.label(RichText::new(info).small().monospace().weak());
        }
        Err(err) => {
            ui.colored_label(Color32::from_rgb(255, 115, 115), RichText::new(err).small());
        }
    }
}

fn mm_badge(ui: &mut egui::Ui, mm: usize) {
    let color = match mm {
        0 => Color32::from_rgb(76, 195, 138),
        1 => Color32::from_rgb(245, 195, 70),
        _ => Color32::from_rgb(255, 130, 90),
    };
    ui.colored_label(color, RichText::new(mm.to_string()).strong());
}

fn render_colored_alignment(
    ui: &mut egui::Ui,
    primer: &str,
    template: &str,
    end_identity: usize,
    is_forward: bool,
) {
    let p_masks = encode_primer(primer);
    let p_bytes = primer.as_bytes();
    let t_bytes = template.as_bytes();
    let len = p_bytes.len().min(t_bytes.len());
    let ei = end_identity.min(len);

    let in_end_id = |i: usize| -> bool {
        if is_forward {
            i >= len - ei
        } else {
            i < ei
        }
    };

    let green = Color32::from_rgb(76, 195, 138);
    let red = Color32::from_rgb(255, 95, 95);
    let normal = ui.visuals().text_color();
    let muted = Color32::GRAY;
    let font = egui::FontId::monospace(13.0);

    let mut primer_job = egui::text::LayoutJob::default();
    let mut mid_job = egui::text::LayoutJob::default();
    let mut tmpl_job = egui::text::LayoutJob::default();

    primer_job.append("Primer:   ", 0.0, egui::TextFormat::simple(font.clone(), muted));
    mid_job.append("          ", 0.0, egui::TextFormat::simple(font.clone(), muted));
    tmpl_job.append("Template: ", 0.0, egui::TextFormat::simple(font.clone(), muted));

    for i in 0..len {
        let matched = bases_match(t_bytes[i], p_masks[i]);
        let color = if !matched {
            red
        } else if in_end_id(i) {
            green
        } else {
            normal
        };
        let mid_ch = if matched { "|" } else { "×" };

        primer_job.append(
            &(p_bytes[i] as char).to_string(),
            0.0,
            egui::TextFormat::simple(font.clone(), color),
        );
        mid_job.append(mid_ch, 0.0, egui::TextFormat::simple(font.clone(), color));
        tmpl_job.append(
            &(t_bytes[i] as char).to_string(),
            0.0,
            egui::TextFormat::simple(font.clone(), color),
        );
    }

    ui.label(primer_job);
    ui.label(mid_job);
    ui.label(tmpl_job);
}

fn format_num(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out.chars().rev().collect()
}
