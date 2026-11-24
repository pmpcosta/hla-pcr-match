//! Downloads and caches official IPD-IMGT/HLA FASTA databases from ANHIG/IMGTHLA on GitHub.

use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::Sender;

const RELEASE_URL: &str =
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/release_version.txt";
const CDS_URL: &str = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_nuc.fasta";
const GEN_ZIP_URL: &str =
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_gen.fasta.zip";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceType {
    Cds,
    Genomic,
}

impl SequenceType {
    pub fn label(self) -> &'static str {
        match self {
            Self::Cds => "CDS database",
            Self::Genomic => "Genomic database",
        }
    }

    pub fn cache_filename(self) -> &'static str {
        match self {
            Self::Cds => "hla_nuc.fasta",
            Self::Genomic => "hla_gen.fasta",
        }
    }
}

#[derive(Debug, Clone)]
pub struct DownloadProgress {
    pub stage: String,
    pub bytes_done: u64,
    pub bytes_total: Option<u64>,
}

#[derive(Debug, Clone, Default)]
pub struct ReleaseInfo {
    pub version: Option<String>,
    pub date: Option<String>,
}

impl ReleaseInfo {
    pub fn display(&self) -> String {
        match (&self.version, &self.date) {
            (Some(v), Some(d)) => format!("{v} ({d})"),
            (Some(v), None) => v.clone(),
            _ => "Local / Unknown release".to_string(),
        }
    }
}

pub fn cache_dir() -> PathBuf {
    if let Some(proj) = directories::ProjectDirs::from("org", "hla-tools", "hla-pcr-match") {
        proj.data_local_dir().to_path_buf()
    } else {
        std::env::current_exe()
            .ok()
            .and_then(|p| p.parent().map(|d| d.join("imgt_hla_cache")))
            .unwrap_or_else(|| PathBuf::from("imgt_hla_cache"))
    }
}

pub fn cached_fasta_path(seq_type: SequenceType) -> PathBuf {
    cache_dir().join(seq_type.cache_filename())
}

pub fn read_cached_release_info() -> ReleaseInfo {
    let path = cache_dir().join("release_version.txt");
    let Ok(text) = fs::read_to_string(path) else {
        return ReleaseInfo::default();
    };
    parse_release_info(&text)
}

pub fn parse_release_info(text: &str) -> ReleaseInfo {
    let mut info = ReleaseInfo::default();
    for line in text.lines() {
        let trimmed = line.trim();
        if let Some(v) = trimmed.strip_prefix("# version:") {
            info.version = Some(v.trim().to_string());
        } else if let Some(d) = trimmed.strip_prefix("# date:") {
            info.date = Some(d.trim().to_string());
        }
    }
    info
}

fn download_file(
    url: &str,
    dest: &Path,
    stage: &str,
    progress_tx: &Sender<DownloadProgress>,
) -> Result<(), String> {
    let resp = ureq::get(url)
        .set("User-Agent", "hla-pcr-match/0.1")
        .call()
        .map_err(|e| format!("HTTP request failed for {url}: {e}"))?;

    let total: Option<u64> = resp
        .header("Content-Length")
        .and_then(|v| v.parse::<u64>().ok());

    let tmp_path = dest.with_extension("part");
    let mut out = File::create(&tmp_path)
        .map_err(|e| format!("Failed to create {}: {e}", tmp_path.display()))?;
    let mut reader = resp.into_reader();
    let mut buf = [0u8; 64 * 1024];
    let mut done: u64 = 0;

    let _ = progress_tx.send(DownloadProgress {
        stage: stage.to_string(),
        bytes_done: 0,
        bytes_total: total,
    });

    loop {
        let n = reader
            .read(&mut buf)
            .map_err(|e| format!("Download read error ({stage}): {e}"))?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])
            .map_err(|e| format!("Disk write error ({stage}): {e}"))?;
        done += n as u64;
        let _ = progress_tx.send(DownloadProgress {
            stage: stage.to_string(),
            bytes_done: done,
            bytes_total: total,
        });
    }
    out.flush()
        .map_err(|e| format!("Flush error ({stage}): {e}"))?;
    drop(out);

    if dest.exists() {
        let _ = fs::remove_file(dest);
    }
    fs::rename(&tmp_path, dest)
        .map_err(|e| format!("Failed to finalize {}: {e}", dest.display()))?;
    Ok(())
}

/// Downloads IMGT/HLA release metadata + requested FASTA(s) into the local cache directory.
pub fn download_imgt_hla(
    download_cds: bool,
    download_genomic: bool,
    progress_tx: &Sender<DownloadProgress>,
) -> Result<ReleaseInfo, String> {
    let dir = cache_dir();
    fs::create_dir_all(&dir)
        .map_err(|e| format!("Failed to create cache dir {}: {e}", dir.display()))?;

    let rel_path = dir.join("release_version.txt");
    download_file(
        RELEASE_URL,
        &rel_path,
        "Fetching IMGT/HLA release metadata...",
        progress_tx,
    )?;
    let info = read_cached_release_info();

    if download_cds {
        let cds_path = cached_fasta_path(SequenceType::Cds);
        download_file(
            CDS_URL,
            &cds_path,
            "Downloading CDS database (hla_nuc.fasta)...",
            progress_tx,
        )?;
    }

    if download_genomic {
        let zip_path = dir.join("hla_gen.fasta.zip");
        download_file(
            GEN_ZIP_URL,
            &zip_path,
            "Downloading Genomic database (hla_gen.fasta.zip)...",
            progress_tx,
        )?;

        let _ = progress_tx.send(DownloadProgress {
            stage: "Extracting hla_gen.fasta from zip...".to_string(),
            bytes_done: 0,
            bytes_total: None,
        });

        let zip_file = File::open(&zip_path)
            .map_err(|e| format!("Failed to open {}: {e}", zip_path.display()))?;
        let mut archive = zip::ZipArchive::new(zip_file)
            .map_err(|e| format!("Failed to read zip archive: {e}"))?;

        let mut extracted = false;
        let gen_dest = cached_fasta_path(SequenceType::Genomic);
        let gen_tmp = gen_dest.with_extension("part");

        for i in 0..archive.len() {
            let mut entry = archive
                .by_index(i)
                .map_err(|e| format!("Zip entry error: {e}"))?;
            let name = entry.name().to_ascii_lowercase();
            if name.ends_with(".fasta") || name.ends_with(".fa") {
                let mut out = File::create(&gen_tmp)
                    .map_err(|e| format!("Failed to create {}: {e}", gen_tmp.display()))?;
                io::copy(&mut entry, &mut out)
                    .map_err(|e| format!("Failed to extract {}: {e}", entry.name()))?;
                out.flush().map_err(|e| format!("Flush error: {e}"))?;
                extracted = true;
                break;
            }
        }
        drop(archive);
        let _ = fs::remove_file(&zip_path);

        if !extracted {
            let _ = fs::remove_file(&gen_tmp);
            return Err("Downloaded hla_gen.fasta.zip did not contain a FASTA file".to_string());
        }
        if gen_dest.exists() {
            let _ = fs::remove_file(&gen_dest);
        }
        fs::rename(&gen_tmp, &gen_dest)
            .map_err(|e| format!("Failed to finalize {}: {e}", gen_dest.display()))?;
    }

    Ok(info)
}
