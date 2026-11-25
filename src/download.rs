//! Downloads and caches official IPD-IMGT/HLA FASTA databases from ANHIG/IMGTHLA on GitHub.

use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::Sender;
use std::time::Duration;

use crate::fasta::is_fasta_zip_entry;

const USER_AGENT: &str = concat!("hla-pcr-match/", env!("CARGO_PKG_VERSION"));
const RELEASE_URL: &str =
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/release_version.txt";
const CDS_URL: &str = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_nuc.fasta";
const GEN_ZIP_URL: &str =
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_gen.fasta.zip";

const CONNECT_TIMEOUT: Duration = Duration::from_secs(30);
const READ_TIMEOUT: Duration = Duration::from_secs(60);

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

    pub fn release_filename(self) -> &'static str {
        match self {
            Self::Cds => "hla_nuc.release.txt",
            Self::Genomic => "hla_gen.release.txt",
        }
    }
}

#[derive(Debug, Clone)]
pub struct DownloadProgress {
    pub stage: String,
    pub bytes_done: u64,
    pub bytes_total: Option<u64>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
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

pub fn read_cached_release_info(seq_type: SequenceType) -> ReleaseInfo {
    let dir = cache_dir();
    let per_type = dir.join(seq_type.release_filename());
    if let Ok(text) = fs::read_to_string(&per_type) {
        return parse_release_info(&text);
    }
    let legacy = dir.join("release_version.txt");
    if let Ok(text) = fs::read_to_string(&legacy) {
        return parse_release_info(&text);
    }
    ReleaseInfo::default()
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

fn build_http_agent() -> ureq::Agent {
    ureq::AgentBuilder::new()
        .timeout_connect(CONNECT_TIMEOUT)
        .timeout_read(READ_TIMEOUT)
        .user_agent(USER_AGENT)
        .build()
}

fn replace_file(tmp_path: &Path, dest: &Path) -> Result<(), String> {
    #[cfg(unix)]
    {
        fs::rename(tmp_path, dest)
            .map_err(|e| format!("Failed to finalize {}: {e}", dest.display()))
    }
    #[cfg(not(unix))]
    {
        if dest.exists() {
            let backup = dest.with_extension("bak");
            if backup.exists() {
                let _ = fs::remove_file(&backup);
            }
            fs::rename(dest, &backup)
                .map_err(|e| format!("Failed to replace {}: {e}", dest.display()))?;
            if let Err(e) = fs::rename(tmp_path, dest) {
                let _ = fs::rename(&backup, dest);
                return Err(format!("Failed to finalize {}: {e}", dest.display()));
            }
            let _ = fs::remove_file(&backup);
            Ok(())
        } else {
            fs::rename(tmp_path, dest)
                .map_err(|e| format!("Failed to finalize {}: {e}", dest.display()))
        }
    }
}

fn write_text_atomically(dest: &Path, text: &str) -> Result<(), String> {
    let tmp_path = dest.with_extension("part");
    let write_res = (|| -> Result<(), String> {
        let mut out = File::create(&tmp_path)
            .map_err(|e| format!("Failed to create {}: {e}", tmp_path.display()))?;
        out.write_all(text.as_bytes())
            .map_err(|e| format!("Disk write error ({}): {e}", dest.display()))?;
        out.flush()
            .map_err(|e| format!("Flush error ({}): {e}", dest.display()))?;
        Ok(())
    })();

    if let Err(e) = write_res {
        let _ = fs::remove_file(&tmp_path);
        return Err(e);
    }

    if let Err(e) = replace_file(&tmp_path, dest) {
        let _ = fs::remove_file(&tmp_path);
        return Err(e);
    }
    Ok(())
}

fn download_file(
    agent: &ureq::Agent,
    url: &str,
    dest: &Path,
    stage: &str,
    progress_tx: &Sender<DownloadProgress>,
    cancel: &AtomicBool,
) -> Result<(), String> {
    if cancel.load(Ordering::Relaxed) {
        return Err("Download cancelled.".to_string());
    }

    let resp = agent
        .get(url)
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

    let stream_res = (|| -> Result<(), String> {
        loop {
            if cancel.load(Ordering::Relaxed) {
                return Err("Download cancelled.".to_string());
            }
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
        Ok(())
    })();

    drop(out);
    if let Err(e) = stream_res {
        let _ = fs::remove_file(&tmp_path);
        return Err(e);
    }

    if let Err(e) = replace_file(&tmp_path, dest) {
        let _ = fs::remove_file(&tmp_path);
        return Err(e);
    }
    Ok(())
}

/// Downloads IMGT/HLA release metadata + requested FASTA(s) into the local cache directory.
pub fn download_imgt_hla(
    download_cds: bool,
    download_genomic: bool,
    progress_tx: &Sender<DownloadProgress>,
    cancel: &AtomicBool,
) -> Result<ReleaseInfo, String> {
    let dir = cache_dir();
    fs::create_dir_all(&dir)
        .map_err(|e| format!("Failed to create cache dir {}: {e}", dir.display()))?;

    let agent = build_http_agent();

    if cancel.load(Ordering::Relaxed) {
        return Err("Download cancelled.".to_string());
    }

    let _ = progress_tx.send(DownloadProgress {
        stage: "Fetching IMGT/HLA release metadata...".to_string(),
        bytes_done: 0,
        bytes_total: None,
    });

    let release_text = agent
        .get(RELEASE_URL)
        .call()
        .map_err(|e| format!("HTTP request failed for {RELEASE_URL}: {e}"))?
        .into_string()
        .map_err(|e| format!("Failed to read IMGT/HLA release metadata: {e}"))?;
    let info = parse_release_info(&release_text);

    if download_cds {
        let cds_path = cached_fasta_path(SequenceType::Cds);
        download_file(
            &agent,
            CDS_URL,
            &cds_path,
            "Downloading CDS database (hla_nuc.fasta)...",
            progress_tx,
            cancel,
        )?;
        write_text_atomically(
            &dir.join(SequenceType::Cds.release_filename()),
            &release_text,
        )?;
    }

    if download_genomic {
        let zip_path = dir.join("hla_gen.fasta.zip");
        download_file(
            &agent,
            GEN_ZIP_URL,
            &zip_path,
            "Downloading Genomic database (hla_gen.fasta.zip)...",
            progress_tx,
            cancel,
        )?;

        let _ = progress_tx.send(DownloadProgress {
            stage: "Extracting hla_gen.fasta from zip...".to_string(),
            bytes_done: 0,
            bytes_total: None,
        });

        let gen_dest = cached_fasta_path(SequenceType::Genomic);
        let gen_tmp = gen_dest.with_extension("part");

        let extract_res = (|| -> Result<(), String> {
            let zip_file = File::open(&zip_path)
                .map_err(|e| format!("Failed to open {}: {e}", zip_path.display()))?;
            let mut archive = zip::ZipArchive::new(zip_file)
                .map_err(|e| format!("Failed to read zip archive: {e}"))?;

            let mut extracted = false;
            let mut buf = [0u8; 64 * 1024];

            for i in 0..archive.len() {
                let mut entry = archive
                    .by_index(i)
                    .map_err(|e| format!("Zip entry error: {e}"))?;
                if is_fasta_zip_entry(entry.name(), entry.is_dir()) {
                    let mut out = File::create(&gen_tmp)
                        .map_err(|e| format!("Failed to create {}: {e}", gen_tmp.display()))?;
                    loop {
                        if cancel.load(Ordering::Relaxed) {
                            return Err("Download cancelled.".to_string());
                        }
                        let n = entry
                            .read(&mut buf)
                            .map_err(|e| format!("Failed to extract {}: {e}", entry.name()))?;
                        if n == 0 {
                            break;
                        }
                        out.write_all(&buf[..n])
                            .map_err(|e| format!("Disk write error during extraction: {e}"))?;
                    }
                    out.flush().map_err(|e| format!("Flush error: {e}"))?;
                    extracted = true;
                    break;
                }
            }

            if !extracted {
                return Err("Downloaded hla_gen.fasta.zip did not contain a FASTA file".to_string());
            }
            Ok(())
        })();

        let _ = fs::remove_file(&zip_path);
        if let Err(e) = extract_res {
            let _ = fs::remove_file(&gen_tmp);
            return Err(e);
        }

        if let Err(e) = replace_file(&gen_tmp, &gen_dest) {
            let _ = fs::remove_file(&gen_tmp);
            return Err(e);
        }
        write_text_atomically(
            &dir.join(SequenceType::Genomic.release_filename()),
            &release_text,
        )?;
    }

    // Keep legacy release_version.txt updated once all requested downloads complete
    let _ = write_text_atomically(&dir.join("release_version.txt"), &release_text);

    Ok(info)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_release_info() {
        let sample = "\
# file: release_version.txt
# date: 2026-04-15
# version: IPD-IMGT/HLA 3.60.0
# origin: https://github.com/ANHIG/IMGTHLA
";
        let info = parse_release_info(sample);
        assert_eq!(info.version.as_deref(), Some("IPD-IMGT/HLA 3.60.0"));
        assert_eq!(info.date.as_deref(), Some("2026-04-15"));
        assert_eq!(info.display(), "IPD-IMGT/HLA 3.60.0 (2026-04-15)");

        let version_only = parse_release_info("# version: 3.60.0\n");
        assert_eq!(version_only.display(), "3.60.0");

        let empty = parse_release_info("");
        assert_eq!(empty.display(), "Local / Unknown release");
    }
}
