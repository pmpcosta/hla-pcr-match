//! Native egui/eframe GUI for HLA Primer Match.

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
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
use crate::fasta::{AlleleDatabase, parse_query_tokens};
use crate::iupac::{
    bases_match, clean_and_validate_primer, encode_primer, primer_degeneracy, reverse_complement,
};
use crate::matcher::{AmpliconHit, MAX_AMPLICONS_PER_ALLELE, MatchParams, MatchResults, run_match};

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

#[derive(Debug, Clone, Copy, Default)]
struct TableScrollState {
    offset_y: f32,
    viewport_h: f32,
    content_h: f32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ScrollKey {
    Up,
    Down,
    PageUp,
    PageDown,
    Home,
    End,
}

struct ThemeColors {
    ok: Color32,
    warn: Color32,
    err: Color32,
    err_bg: Color32,
    info: Color32,
    mm_zero: Color32,
    mm_one: Color32,
    mm_multi: Color32,
    /// (background, foreground) for mismatched bases in a site with 1 mismatch
    site_mm1: (Color32, Color32),
    /// (background, foreground) for mismatched bases in a site with 2 mismatches
    site_mm2: (Color32, Color32),
    /// (background, foreground) for mismatched bases in a site with 3+ mismatches
    site_mm3_plus: (Color32, Color32),
    muted: Color32,
}

fn theme_colors(dark_mode: bool) -> ThemeColors {
    if dark_mode {
        ThemeColors {
            ok: Color32::from_rgb(76, 195, 138),
            warn: Color32::from_rgb(245, 179, 66),
            err: Color32::from_rgb(255, 115, 115),
            err_bg: Color32::from_rgba_unmultiplied(200, 50, 50, 35),
            info: Color32::from_rgb(130, 180, 255),
            mm_zero: Color32::from_rgb(76, 195, 138),
            mm_one: Color32::from_rgb(245, 195, 70),
            mm_multi: Color32::from_rgb(255, 130, 90),
            site_mm1: (
                Color32::from_rgb(38, 96, 62),
                Color32::from_rgb(220, 252, 231),
            ),
            site_mm2: (
                Color32::from_rgb(120, 92, 22),
                Color32::from_rgb(254, 249, 195),
            ),
            site_mm3_plus: (
                Color32::from_rgb(138, 42, 42),
                Color32::from_rgb(254, 226, 226),
            ),
            muted: Color32::GRAY,
        }
    } else {
        ThemeColors {
            ok: Color32::from_rgb(18, 128, 72),
            warn: Color32::from_rgb(168, 100, 0),
            err: Color32::from_rgb(190, 34, 34),
            err_bg: Color32::from_rgba_unmultiplied(210, 40, 40, 24),
            info: Color32::from_rgb(22, 92, 185),
            mm_zero: Color32::from_rgb(18, 128, 72),
            mm_one: Color32::from_rgb(168, 100, 0),
            mm_multi: Color32::from_rgb(196, 56, 20),
            site_mm1: (
                Color32::from_rgb(187, 247, 208),
                Color32::from_rgb(20, 83, 45),
            ),
            site_mm2: (
                Color32::from_rgb(254, 240, 138),
                Color32::from_rgb(113, 63, 18),
            ),
            site_mm3_plus: (
                Color32::from_rgb(254, 202, 202),
                Color32::from_rgb(127, 29, 29),
            ),
            muted: Color32::from_rgb(100, 100, 100),
        }
    }
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
    cds_release: ReleaseInfo,
    genomic_release: ReleaseInfo,
    cds_db: Option<Arc<AlleleDatabase>>,
    genomic_db: Option<Arc<AlleleDatabase>>,
    custom_db: Option<Arc<AlleleDatabase>>,
    loading_db: bool,

    // Background worker state (persistent channel avoids losing in-flight events)
    bg_tx: Sender<BgEvent>,
    bg_rx: Receiver<BgEvent>,
    download_progress: Option<DownloadProgress>,
    download_cancel: Option<Arc<AtomicBool>>,
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
    filtered_coverage_indices: Vec<usize>,
    filtered_miss_indices: Vec<usize>,
    indices_dirty: bool,
    matches_scroll: TableScrollState,
    coverage_scroll: TableScrollState,
    missed_scroll: TableScrollState,

    // Feedback banners
    error_msg: Option<String>,
    status_msg: Option<String>,
    dark_mode: bool,
}

impl HlaPrimerApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let dark_mode = true;
        cc.egui_ctx.set_visuals(egui::Visuals::dark());

        let (bg_tx, bg_rx) = mpsc::channel();
        let cds_release = read_cached_release_info(SequenceType::Cds);
        let genomic_release = read_cached_release_info(SequenceType::Genomic);

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
            cds_release,
            genomic_release,
            cds_db: None,
            genomic_db: None,
            custom_db: None,
            loading_db: false,

            bg_tx,
            bg_rx,
            download_progress: None,
            download_cancel: None,
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
            filtered_coverage_indices: Vec::new(),
            filtered_miss_indices: Vec::new(),
            indices_dirty: false,
            matches_scroll: TableScrollState::default(),
            coverage_scroll: TableScrollState::default(),
            missed_scroll: TableScrollState::default(),

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

    fn active_release_info(&self) -> &ReleaseInfo {
        match self.seq_type {
            SequenceType::Cds => &self.cds_release,
            SequenceType::Genomic => &self.genomic_release,
        }
    }

    fn is_busy(&self) -> bool {
        self.loading_db || self.download_progress.is_some() || self.matching_in_progress
    }

    fn start_load_cached_db(&mut self, seq_type: SequenceType, ctx: egui::Context) {
        if self.is_busy() {
            return;
        }
        let path = cached_fasta_path(seq_type);
        if !path.exists() {
            return;
        }
        self.loading_db = true;
        self.error_msg = None;
        self.status_msg = Some(format!("Loading {}...", seq_type.label()));

        let tx = self.bg_tx.clone();
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
        if self.is_busy() {
            return;
        }
        self.loading_db = true;
        self.error_msg = None;
        let fname = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("custom.fasta")
            .to_string();
        self.status_msg = Some(format!("Loading {fname}..."));

        let tx = self.bg_tx.clone();

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
        let cancel = Arc::new(AtomicBool::new(false));
        self.error_msg = None;
        self.status_msg = None;
        self.download_cancel = Some(Arc::clone(&cancel));
        self.download_progress = Some(DownloadProgress {
            stage: "Connecting to GitHub (ANHIG/IMGTHLA)...".to_string(),
            bytes_done: 0,
            bytes_total: None,
        });

        let tx = self.bg_tx.clone();

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

            let result = download_imgt_hla(cds, genomic, &prog_tx, &cancel);
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

        // Release focus from any TextEdit so arrow / PageUp / PageDown keys scroll results directly.
        ctx.memory_mut(|m| {
            if let Some(id) = m.focused() {
                m.surrender_focus(id);
            }
        });

        let tx = self.bg_tx.clone();

        thread::spawn(move || {
            let res = run_match(db, &params);
            let _ = tx.send(BgEvent::MatchFinished(res));
            ctx.request_repaint();
        });
    }

    fn poll_background(&mut self, ctx: &egui::Context) {
        let mut events = Vec::new();
        while let Ok(ev) = self.bg_rx.try_recv() {
            events.push(ev);
        }

        for ev in events {
            match ev {
                BgEvent::DownloadProgress(p) => {
                    if self.download_progress.is_some() {
                        self.download_progress = Some(p);
                    }
                }
                BgEvent::DownloadFinished(res, downloaded_cds, downloaded_gen) => {
                    self.download_progress = None;
                    self.download_cancel = None;
                    // Refresh per-DB cached release info even on partial success
                    self.cds_release = read_cached_release_info(SequenceType::Cds);
                    self.genomic_release = read_cached_release_info(SequenceType::Genomic);
                    match res {
                        Ok(_info) => {
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
                            if err == "Download cancelled." {
                                self.status_msg = Some(err);
                            } else {
                                self.error_msg = Some(format!("Download failed: {err}"));
                            }
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
            self.filtered_coverage_indices.clear();
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
                    || contains_ascii_ci(&h.allele, &needle)
                    || contains_ascii_ci(&h.locus, &needle)
                    || contains_ascii_ci(&h.accession, &needle)
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
                HitSortCol::TotalMm => (ha.fwd_mismatches + ha.rev_mismatches)
                    .cmp(&(hb.fwd_mismatches + hb.rev_mismatches)),
            };
            let ordered = if asc { ord } else { ord.reverse() };
            ordered.then_with(|| ha.allele.cmp(&hb.allele))
        });

        self.filtered_coverage_indices = res
            .locus_coverage
            .iter()
            .enumerate()
            .filter(|(_, c)| needle.is_empty() || contains_ascii_ci(&c.locus, &needle))
            .map(|(i, _)| i)
            .collect();

        self.filtered_miss_indices = res
            .missed
            .iter()
            .enumerate()
            .filter(|(_, m)| {
                needle.is_empty()
                    || contains_ascii_ci(&m.allele, &needle)
                    || contains_ascii_ci(&m.locus, &needle)
                    || contains_ascii_ci(&m.reason.label(), &needle)
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
                MissSortCol::Reason => ma.reason.sort_key().cmp(&mb.reason.sort_key()),
                MissSortCol::FwdHits => ma.fwd_hits.cmp(&mb.fwd_hits),
                MissSortCol::RevHits => ma.rev_hits.cmp(&mb.rev_hits),
            };
            let ordered = if masc { ord } else { ord.reverse() };
            ordered.then_with(|| ma.allele.cmp(&mb.allele))
        });

        self.indices_dirty = false;
    }

    fn load_example(&mut self, ctx: &egui::Context) {
        self.query = "HLA-A".to_string();
        // Real conserved HLA-A exon 1 -> exon 2 CDS primers from A*01:01:01:01
        self.primer_fwd = "ATGGCCGTCATGGCGCCCCG".to_string();
        self.primer_rvs = "CTGGTCCCAATACTCCGG".to_string();
        self.end_identity = 3;
        self.max_mismatch = 3;
        self.flip_reverse = true;
        self.product_size_min = -1;
        self.product_size_max = -1;
        self.error_msg = None;

        if !self.is_busy() {
            self.source_mode = DataSourceMode::ImgtHla;
            self.seq_type = SequenceType::Cds;
            if self.cds_db.is_none() && cached_fasta_path(SequenceType::Cds).exists() {
                self.start_load_cached_db(SequenceType::Cds, ctx.clone());
            }
        }
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
        let colors = theme_colors(self.dark_mode);
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
                            self.active_release_info().display()
                        ),
                        DataSourceMode::CustomFasta => {
                            format!(
                                "Custom: {} • {} alleles",
                                db.source_label,
                                format_num(db.len())
                            )
                        }
                    };
                    ui.colored_label(colors.ok, badge);
                } else {
                    ui.colored_label(colors.warn, format!("{} not loaded", self.seq_type.label()));
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
                        if ui.button("Open FASTA...").clicked()
                            && let Some(path) = rfd::FileDialog::new()
                                .add_filter("FASTA / ZIP", &["fasta", "fa", "fna", "zip"])
                                .pick_file()
                        {
                            self.start_load_custom_fasta(path, ctx.clone());
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
                let frac = prog.bytes_total.map(|t| {
                    if t > 0 {
                        prog.bytes_done as f32 / t as f32
                    } else {
                        0.0
                    }
                });
                let text = match prog.bytes_total {
                    Some(total) => format!(
                        "{} ({:.1} / {:.1} MB)",
                        prog.stage,
                        prog.bytes_done as f64 / 1_048_576.0,
                        total as f64 / 1_048_576.0
                    ),
                    None => prog.stage.clone(),
                };
                let mut cancel_clicked = false;
                ui.horizontal(|ui| {
                    if let Some(f) = frac {
                        ui.add(
                            egui::ProgressBar::new(f)
                                .desired_width((ui.available_width() - 74.0).max(120.0))
                                .text(text)
                                .animate(true),
                        );
                    } else {
                        ui.spinner();
                        ui.label(text);
                    }
                    if ui.small_button("Cancel").clicked() {
                        cancel_clicked = true;
                    }
                });
                if cancel_clicked && let Some(ref flag) = self.download_cancel {
                    flag.store(true, Ordering::Relaxed);
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
        let colors = theme_colors(self.dark_mode);
        let busy = self.is_busy();

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
                            "Choose specific loci, alleles, or accessions separated by a comma or newline.\n\
                             Examples: HLA-A, DRB1, B*07:02, A*01, MICA, HLA:HLA00001.\n\
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
                                let already_present = parse_query_tokens(&self.query)
                                    .iter()
                                    .any(|t| t.raw.eq_ignore_ascii_case(loc));
                                if !already_present {
                                    if self.query.trim().is_empty() {
                                        self.query = loc.to_string();
                                    } else {
                                        self.query.push_str(", ");
                                        self.query.push_str(loc);
                                    }
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
                            "5'→3' sequence of the forward primer. IUPAC degenerate bases (R, Y, S, W, K, M, B, D, H, V, N) are supported; optional 5'- / -3' wrappers and spaces are stripped automatically.",
                        );
                    });
                    ui.add(
                        egui::TextEdit::singleline(&mut self.primer_fwd)
                            .desired_width(f32::INFINITY)
                            .font(TextStyle::Monospace)
                            .hint_text("5'-ATGGCCGTCATGGCGCCCCG-3'"),
                    );
                    render_primer_status(ui, &self.primer_fwd, false, &colors);

                    ui.add_space(8.0);

                    // 3. Reverse primer
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Reverse primer:").strong());
                        help_icon(
                            ui,
                            "Sequence of the reverse primer (IUPAC degenerate bases supported; optional 5'- / -3' wrappers are stripped).\n\
                             When 'Flip reverse' is Yes (default), enter it in standard 5'→3' oligo order and it will be reverse-complemented automatically.",
                        );
                    });
                    ui.add(
                        egui::TextEdit::singleline(&mut self.primer_rvs)
                            .desired_width(f32::INFINITY)
                            .font(TextStyle::Monospace)
                            .hint_text("5'-CTGGTCCCAATACTCCGG-3'"),
                    );
                    render_primer_status(ui, &self.primer_rvs, self.flip_reverse, &colors);

                    ui.add_space(10.0);

                    // 4. Type (CDS database / Genomic database)
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Type:").strong());
                        help_icon(
                            ui,
                            "Select coding DNA (CDS: hla_nuc.fasta) or genomic DNA (Genomic: hla_gen.fasta including introns/UTRs).",
                        );
                        let prev_type = self.seq_type;
                        ui.add_enabled_ui(!busy, |ui| {
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
                            ui.colored_label(colors.info, "Using custom FASTA file.");
                            if ui
                                .add_enabled(!busy, egui::Button::new("Switch to IMGT/HLA").small())
                                .clicked()
                            {
                                self.source_mode = DataSourceMode::ImgtHla;
                                if self.active_db().is_none()
                                    && cached_fasta_path(self.seq_type).exists()
                                {
                                    self.start_load_cached_db(self.seq_type, ctx.clone());
                                }
                            }
                        });
                    } else if self.active_db().is_none() && !busy {
                        ui.horizontal(|ui| {
                            ui.colored_label(
                                colors.warn,
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

                    // Action buttons (Ctrl+Enter on Linux/Windows, Cmd+Enter on macOS)
                    let cmd_enter =
                        ctx.input(|i| i.modifiers.command && i.key_pressed(egui::Key::Enter));
                    ui.horizontal(|ui| {
                        let submit_btn = egui::Button::new(
                            RichText::new("  Submit  ").strong().color(Color32::WHITE),
                        )
                        .fill(Color32::from_rgb(38, 128, 214))
                        .min_size(egui::vec2(110.0, 32.0));

                        let clicked = ui
                            .add_enabled(!busy, submit_btn)
                            .on_hover_text("Run primer match (Ctrl+Enter / Cmd+Enter)")
                            .clicked();
                        if clicked || (cmd_enter && !busy) {
                            self.start_match(ctx.clone());
                        }

                        if ui
                            .add_sized([76.0, 32.0], egui::Button::new("Example"))
                            .clicked()
                        {
                            self.load_example(ctx);
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
                            .fill(colors.err_bg)
                            .show(ui, |ui| {
                                ui.colored_label(colors.err, err);
                            });
                    }
                });
            });
    }

    fn render_results_panel(&mut self, ctx: &egui::Context) {
        let colors = theme_colors(self.dark_mode);

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
                        if ui.button("Or Open a Local FASTA / ZIP File...").clicked()
                            && let Some(path) = rfd::FileDialog::new()
                                .add_filter("FASTA / ZIP", &["fasta", "fa", "fna", "zip"])
                                .pick_file()
                        {
                            self.start_load_custom_fasta(path, ctx.clone());
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
            let (
                matched_u,
                total_s,
                hit_cnt,
                miss_cnt,
                locus_cnt,
                elapsed,
                unmatched_q,
                truncated_hits,
                source_label,
            ) = {
                let res = self.results.as_ref().unwrap();
                (
                    res.matched_unique_alleles,
                    res.total_searched,
                    res.hits.len(),
                    res.missed.len(),
                    res.locus_coverage.len(),
                    res.elapsed_ms,
                    res.unmatched_queries.clone(),
                    res.truncated_hits,
                    res.db.source_label.clone(),
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
                            colors.ok
                        } else if cov_pct > 0.0 {
                            colors.warn
                        } else {
                            colors.err
                        }),
                    );
                    ui.separator();
                    ui.label(format!("Amplicons: {}", format_num(hit_cnt)));
                    ui.separator();
                    ui.label(format!("Missed: {}", format_num(miss_cnt)));
                    ui.separator();
                    ui.label(RichText::new(format!("{source_label} • {elapsed} ms")).weak());

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
                        colors.warn,
                        format!(
                            "Warning — unrecognized query term(s): {}",
                            unmatched_q.join(", ")
                        ),
                    );
                }

                if truncated_hits {
                    ui.add_space(4.0);
                    ui.colored_label(
                        colors.warn,
                        format!(
                            "Warning — one or more alleles exceeded the per-allele limit ({MAX_AMPLICONS_PER_ALLELE} amplicons) and were truncated. Tighten primer specificity or Product Size bounds."
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
                            .id(egui::Id::new("results_filter_input"))
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
                (ui.available_height() - 195.0).max(160.0)
            } else {
                ui.available_height().max(160.0)
            };

            ui.allocate_ui_with_layout(
                egui::vec2(ui.available_width(), table_height),
                egui::Layout::top_down(egui::Align::Min),
                |ui| match self.active_tab {
                    ActiveTab::Matches => self.render_matches_table(ui, &colors),
                    ActiveTab::Coverage => self.render_coverage_table(ui),
                    ActiveTab::Missed => self.render_missed_table(ui),
                },
            );

            if show_inspector {
                ui.separator();
                self.render_alignment_inspector(ui, &colors);
            }
        });
    }

    fn render_matches_table(&mut self, ui: &mut egui::Ui, colors: &ThemeColors) {
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
        let fwd_masks = encode_primer(&res.fwd_primer_used);
        let rev_masks = encode_primer(&res.rev_motif_used);

        const ROW_HEIGHT: f32 = 20.0;
        let mut scroll_to_row: Option<usize> = None;
        let mut scroll_offset: Option<f32> = None;

        if !ui.ctx().wants_keyboard_input() && ui.ctx().input(|i| i.key_pressed(egui::Key::Escape))
        {
            self.selected_hit_idx = None;
        }

        if let Some(key) = read_scroll_key(ui.ctx()) {
            let row_count = self.filtered_hit_indices.len();
            let cur_selected_row = self
                .selected_hit_idx
                .and_then(|sel| self.filtered_hit_indices.iter().position(|&idx| idx == sel));

            if let Some(cur_row) = cur_selected_row {
                let page_rows = if self.matches_scroll.viewport_h > 0.0 {
                    ((self.matches_scroll.viewport_h / ROW_HEIGHT).floor() as usize)
                        .saturating_sub(1)
                        .max(1)
                } else {
                    10
                };
                let next_row = match key {
                    ScrollKey::Down => (cur_row + 1).min(row_count - 1),
                    ScrollKey::Up => cur_row.saturating_sub(1),
                    ScrollKey::PageDown => (cur_row + page_rows).min(row_count - 1),
                    ScrollKey::PageUp => cur_row.saturating_sub(page_rows),
                    ScrollKey::Home => 0,
                    ScrollKey::End => row_count - 1,
                };
                self.selected_hit_idx = Some(self.filtered_hit_indices[next_row]);
                scroll_to_row = Some(next_row);
            } else {
                scroll_offset = Some(next_scroll_offset(self.matches_scroll, ROW_HEIGHT, key));
            }
        }

        let mut builder = TableBuilder::new(ui)
            .id_salt("matches_table")
            .striped(true)
            .resizable(true)
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
            .sense(egui::Sense::click())
            .column(Column::initial(155.0).at_least(100.0)) // Allele
            .column(Column::initial(65.0).at_least(50.0)) // Start
            .column(Column::initial(65.0).at_least(50.0)) // End
            .column(Column::initial(70.0).at_least(55.0)) // Length
            .column(Column::initial(75.0).at_least(55.0)) // Total MM
            .column(Column::remainder().at_least(180.0)); // Binding sites preview

        if let Some(row) = scroll_to_row {
            builder = builder.scroll_to_row(row, None);
        } else if let Some(offset) = scroll_offset {
            builder = builder.vertical_scroll_offset(offset);
        }

        let scroll_out = builder
            .header(22.0, |mut header| {
                let cols = [
                    ("Allele", HitSortCol::Allele, None),
                    (
                        "Start",
                        HitSortCol::Start,
                        Some("0-based inclusive start coordinate"),
                    ),
                    (
                        "End",
                        HitSortCol::End,
                        Some("0-based exclusive end coordinate"),
                    ),
                    ("Length", HitSortCol::Length, Some("Amplicon length in bp")),
                    ("Total MM", HitSortCol::TotalMm, None),
                ];
                for (title, col, tip) in cols {
                    header.col(|ui| {
                        let arrow = if self.hit_sort_col == col {
                            if self.hit_sort_asc { " ▲" } else { " ▼" }
                        } else {
                            ""
                        };
                        let mut resp = ui.button(RichText::new(format!("{title}{arrow}")).strong());
                        if let Some(t) = tip {
                            resp = resp.on_hover_text(t);
                        }
                        if resp.clicked() {
                            sort_change = Some(col);
                        }
                    });
                }
                header.col(|ui| {
                    ui.label(RichText::new("Fwd / Rev Template Sites").strong())
                        .on_hover_text(
                            "Mismatched base background color by per-site mismatch count:\n\
                             Green = 1 mismatch • Yellow = 2 mismatches • Red = 3+ mismatches",
                        );
                });
            })
            .body(|body| {
                let row_count = self.filtered_hit_indices.len();
                body.rows(ROW_HEIGHT, row_count, |mut row| {
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
                        mm_badge(ui, hit.fwd_mismatches + hit.rev_mismatches, colors);
                    });
                    row.col(|ui| {
                        render_template_sites_cell(ui, hit, &fwd_masks, &rev_masks, colors);
                    });

                    if row.response().clicked() {
                        newly_selected = Some(hit_idx);
                    }
                });
            });

        self.matches_scroll = TableScrollState {
            offset_y: scroll_out.state.offset.y,
            viewport_h: scroll_out.inner_rect.height(),
            content_h: scroll_out.content_size.y,
        };

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

    fn render_coverage_table(&mut self, ui: &mut egui::Ui) {
        let Some(ref res) = self.results else {
            return;
        };

        const ROW_HEIGHT: f32 = 22.0;
        let scroll_offset = read_scroll_key(ui.ctx())
            .map(|key| next_scroll_offset(self.coverage_scroll, ROW_HEIGHT, key));

        let mut builder = TableBuilder::new(ui)
            .id_salt("coverage_table")
            .striped(true)
            .resizable(true)
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
            .column(Column::initial(110.0))
            .column(Column::initial(85.0))
            .column(Column::initial(85.0))
            .column(Column::initial(85.0))
            .column(Column::remainder().at_least(180.0));

        if let Some(offset) = scroll_offset {
            builder = builder.vertical_scroll_offset(offset);
        }

        let scroll_out = builder
            .header(22.0, |mut header| {
                for title in ["Locus", "Matched", "Missed", "Total", "Coverage"] {
                    header.col(|ui| {
                        ui.label(RichText::new(title).strong());
                    });
                }
            })
            .body(|body| {
                let count = self.filtered_coverage_indices.len();
                body.rows(ROW_HEIGHT, count, |mut row| {
                    let c = &res.locus_coverage[self.filtered_coverage_indices[row.index()]];
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

        self.coverage_scroll = TableScrollState {
            offset_y: scroll_out.state.offset.y,
            viewport_h: scroll_out.inner_rect.height(),
            content_h: scroll_out.content_size.y,
        };
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

        const ROW_HEIGHT: f32 = 20.0;
        let scroll_offset = read_scroll_key(ui.ctx())
            .map(|key| next_scroll_offset(self.missed_scroll, ROW_HEIGHT, key));

        let mut sort_change: Option<MissSortCol> = None;

        let mut builder = TableBuilder::new(ui)
            .id_salt("missed_table")
            .striped(true)
            .resizable(true)
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
            .column(Column::initial(160.0))
            .column(Column::initial(80.0))
            .column(Column::remainder().at_least(220.0))
            .column(Column::initial(80.0))
            .column(Column::initial(80.0));

        if let Some(offset) = scroll_offset {
            builder = builder.vertical_scroll_offset(offset);
        }

        let scroll_out = builder
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
                body.rows(ROW_HEIGHT, count, |mut row| {
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

        self.missed_scroll = TableScrollState {
            offset_y: scroll_out.state.offset.y,
            viewport_h: scroll_out.inner_rect.height(),
            content_h: scroll_out.content_size.y,
        };

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

    fn render_alignment_inspector(&mut self, ui: &mut egui::Ui, colors: &ThemeColors) {
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
                if ui.small_button("Copy Amplicon DNA").clicked()
                    && let Some(rec) = res.db.records.get(hit.allele_index)
                    && hit.end <= rec.sequence.len()
                    && hit.start < hit.end
                {
                    copy_seq = Some(
                        String::from_utf8_lossy(&rec.sequence[hit.start..hit.end]).into_owned(),
                    );
                }
            });
        });

        ui.add_space(4.0);
        egui::ScrollArea::horizontal()
            .id_salt("alignment_inspector_scroll")
            .show(ui, |ui| {
                ui.vertical(|ui| {
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
                                colors,
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
                                colors,
                            );
                        });
                    });
                    // Keep space below alignment frames so the horizontal scrollbar never overlaps text.
                    ui.add_space(18.0);
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
        let Some(ref res) = self.results else { return };
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("hla_pcr_amplicons.fasta")
            .add_filter("FASTA", &["fasta", "fa"])
            .save_file()
        {
            match export_amplicons_fasta(&path, res) {
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

        // TopBottomPanel must be added before SidePanel / CentralPanel in egui.
        if let Some(ref msg) = self.status_msg {
            egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(msg).small().weak());
                });
            });
        }

        self.render_form_panel(ctx);
        self.render_results_panel(ctx);
    }
}

fn help_icon(ui: &mut egui::Ui, tooltip: &str) {
    ui.label(RichText::new("(?)").small().weak())
        .on_hover_text(tooltip);
}

fn render_primer_status(ui: &mut egui::Ui, raw: &str, show_rc: bool, colors: &ThemeColors) {
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
            ui.colored_label(colors.err, RichText::new(err).small());
        }
    }
}

fn mm_badge(ui: &mut egui::Ui, mm: usize, colors: &ThemeColors) {
    let color = match mm {
        0 => colors.mm_zero,
        1 => colors.mm_one,
        _ => colors.mm_multi,
    };
    ui.colored_label(color, RichText::new(mm.to_string()).strong());
}

fn append_highlighted_site(
    job: &mut egui::text::LayoutJob,
    template: &str,
    primer_masks: &[u8],
    site_mismatches: usize,
    font: &egui::FontId,
    normal_color: Color32,
    colors: &ThemeColors,
) {
    let normal_fmt = egui::TextFormat::simple(font.clone(), normal_color);
    if site_mismatches == 0 || primer_masks.is_empty() {
        job.append(template, 0.0, normal_fmt);
        return;
    }

    let (mm_bg, mm_fg) = match site_mismatches {
        1 => colors.site_mm1,
        2 => colors.site_mm2,
        _ => colors.site_mm3_plus,
    };
    let mm_fmt = egui::TextFormat {
        font_id: font.clone(),
        color: mm_fg,
        background: mm_bg,
        ..Default::default()
    };

    let mut run_start = 0usize;
    let mut run_is_mm: Option<bool> = None;

    for (pos, (byte_idx, ch)) in template.char_indices().enumerate() {
        let is_mm = if pos < primer_masks.len() {
            !(ch.is_ascii() && bases_match(ch as u8, primer_masks[pos]))
        } else {
            false
        };

        match run_is_mm {
            None => run_is_mm = Some(is_mm),
            Some(prev) if prev != is_mm => {
                let fmt = if prev {
                    mm_fmt.clone()
                } else {
                    normal_fmt.clone()
                };
                job.append(&template[run_start..byte_idx], 0.0, fmt);
                run_start = byte_idx;
                run_is_mm = Some(is_mm);
            }
            Some(_) => {}
        }
    }

    if let Some(prev) = run_is_mm
        && run_start < template.len()
    {
        let fmt = if prev { mm_fmt } else { normal_fmt };
        job.append(&template[run_start..], 0.0, fmt);
    }
}

fn render_template_sites_cell(
    ui: &mut egui::Ui,
    hit: &AmpliconHit,
    fwd_masks: &[u8],
    rev_masks: &[u8],
    colors: &ThemeColors,
) {
    let font = TextStyle::Monospace.resolve(ui.style());
    let normal_color = ui.visuals().text_color();
    let mut job = egui::text::LayoutJob::default();

    append_highlighted_site(
        &mut job,
        &hit.fwd_template,
        fwd_masks,
        hit.fwd_mismatches,
        &font,
        normal_color,
        colors,
    );
    job.append(
        " … ",
        0.0,
        egui::TextFormat::simple(font.clone(), colors.muted),
    );
    append_highlighted_site(
        &mut job,
        &hit.rev_template,
        rev_masks,
        hit.rev_mismatches,
        &font,
        normal_color,
        colors,
    );

    ui.label(job);
}

fn render_colored_alignment(
    ui: &mut egui::Ui,
    primer: &str,
    template: &str,
    end_identity: usize,
    is_forward: bool,
    colors: &ThemeColors,
) {
    let p_masks = encode_primer(primer);
    let p_bytes = primer.as_bytes();
    let t_bytes = template.as_bytes();
    let len = p_bytes.len().min(t_bytes.len());
    let ei = end_identity.min(len);

    let in_end_id = |i: usize| -> bool { if is_forward { i >= len - ei } else { i < ei } };

    let green = colors.ok;
    let red = colors.err;
    let normal = ui.visuals().text_color();
    let muted = colors.muted;
    let font = egui::FontId::monospace(13.0);

    let mut primer_job = egui::text::LayoutJob::default();
    let mut mid_job = egui::text::LayoutJob::default();
    let mut tmpl_job = egui::text::LayoutJob::default();

    primer_job.append(
        "Primer:   ",
        0.0,
        egui::TextFormat::simple(font.clone(), muted),
    );
    mid_job.append(
        "          ",
        0.0,
        egui::TextFormat::simple(font.clone(), muted),
    );
    tmpl_job.append(
        "Template: ",
        0.0,
        egui::TextFormat::simple(font.clone(), muted),
    );

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

fn read_scroll_key(ctx: &egui::Context) -> Option<ScrollKey> {
    if ctx.wants_keyboard_input() {
        return None;
    }
    ctx.input(|i| {
        if i.key_pressed(egui::Key::ArrowDown) {
            Some(ScrollKey::Down)
        } else if i.key_pressed(egui::Key::ArrowUp) {
            Some(ScrollKey::Up)
        } else if i.key_pressed(egui::Key::PageDown) {
            Some(ScrollKey::PageDown)
        } else if i.key_pressed(egui::Key::PageUp) {
            Some(ScrollKey::PageUp)
        } else if i.key_pressed(egui::Key::Home) {
            Some(ScrollKey::Home)
        } else if i.key_pressed(egui::Key::End) {
            Some(ScrollKey::End)
        } else {
            None
        }
    })
}

fn next_scroll_offset(state: TableScrollState, row_height: f32, key: ScrollKey) -> f32 {
    let max_offset = (state.content_h - state.viewport_h).max(0.0);
    let step = row_height * 3.0;
    let page = (state.viewport_h - row_height).max(row_height);
    let target = match key {
        ScrollKey::Down => state.offset_y + step,
        ScrollKey::Up => state.offset_y - step,
        ScrollKey::PageDown => state.offset_y + page,
        ScrollKey::PageUp => state.offset_y - page,
        ScrollKey::Home => 0.0,
        ScrollKey::End => max_offset,
    };
    target.clamp(0.0, max_offset)
}

fn contains_ascii_ci(haystack: &str, needle_lower: &str) -> bool {
    if needle_lower.is_empty() {
        return true;
    }
    let h = haystack.as_bytes();
    let n = needle_lower.as_bytes();
    if n.len() > h.len() {
        return false;
    }
    h.windows(n.len()).any(|w| {
        w.iter()
            .zip(n.iter())
            .all(|(&hb, &nb)| hb.to_ascii_lowercase() == nb)
    })
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
