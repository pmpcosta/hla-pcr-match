//! CSV, TSV, coverage summary, and amplicon FASTA export utilities.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::fasta::AlleleDatabase;
use crate::matcher::MatchResults;

fn escape_field(value: &str, sep: char) -> String {
    if value.contains(sep) || value.contains('"') || value.contains('\n') || value.contains('\r') {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

pub fn export_hits_delimited(
    path: &Path,
    results: &MatchResults,
    sep: char,
) -> Result<usize, String> {
    let file =
        File::create(path).map_err(|e| format!("Failed to create {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "Allele{s}Accession{s}Locus{s}Start{s}End{s}Length{s}FwdMismatches{s}RevMismatches{s}TotalMismatches{s}FwdTemplate{s}RevTemplate",
        s = sep
    )
    .map_err(|e| format!("Write error: {e}"))?;

    for h in &results.hits {
        writeln!(
            w,
            "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
            escape_field(&h.allele, sep),
            escape_field(&h.accession, sep),
            escape_field(&h.locus, sep),
            h.start,
            h.end,
            h.length,
            h.fwd_mismatches,
            h.rev_mismatches,
            h.fwd_mismatches + h.rev_mismatches,
            escape_field(&h.fwd_template, sep),
            escape_field(&h.rev_template, sep),
            s = sep
        )
        .map_err(|e| format!("Write error: {e}"))?;
    }

    w.flush().map_err(|e| format!("Flush error: {e}"))?;
    Ok(results.hits.len())
}

pub fn export_missed_csv(path: &Path, results: &MatchResults) -> Result<usize, String> {
    let file =
        File::create(path).map_err(|e| format!("Failed to create {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "Allele,Accession,Locus,Reason,FwdBindingSites,RevBindingSites")
        .map_err(|e| format!("Write error: {e}"))?;

    for m in &results.missed {
        writeln!(
            w,
            "{},{},{},{},{},{}",
            escape_field(&m.allele, ','),
            escape_field(&m.accession, ','),
            escape_field(&m.locus, ','),
            escape_field(&m.reason.label(), ','),
            m.fwd_hits,
            m.rev_hits,
        )
        .map_err(|e| format!("Write error: {e}"))?;
    }

    w.flush().map_err(|e| format!("Flush error: {e}"))?;
    Ok(results.missed.len())
}

pub fn export_coverage_csv(path: &Path, results: &MatchResults) -> Result<usize, String> {
    let file =
        File::create(path).map_err(|e| format!("Failed to create {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "Locus,MatchedAlleles,MissedAlleles,TotalAlleles,CoveragePct"
    )
    .map_err(|e| format!("Write error: {e}"))?;

    for c in &results.locus_coverage {
        writeln!(
            w,
            "{},{},{},{},{:.2}",
            escape_field(&c.locus, ','),
            c.matched_alleles,
            c.missed_alleles,
            c.total_alleles,
            c.coverage_pct,
        )
        .map_err(|e| format!("Write error: {e}"))?;
    }

    w.flush().map_err(|e| format!("Flush error: {e}"))?;
    Ok(results.locus_coverage.len())
}

pub fn export_amplicons_fasta(
    path: &Path,
    db: &AlleleDatabase,
    results: &MatchResults,
) -> Result<usize, String> {
    let file =
        File::create(path).map_err(|e| format!("Failed to create {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);

    let mut count = 0usize;
    for h in &results.hits {
        if let Some(rec) = db.records.get(h.allele_index) {
            if h.end <= rec.sequence.len() && h.start < h.end {
                let amplicon = &rec.sequence[h.start..h.end];
                writeln!(
                    w,
                    ">{} {} start={} end={} len={} fwd_mm={} rev_mm={}",
                    h.allele,
                    h.accession,
                    h.start,
                    h.end,
                    h.length,
                    h.fwd_mismatches,
                    h.rev_mismatches
                )
                .map_err(|e| format!("Write error: {e}"))?;

                for chunk in amplicon.chunks(80) {
                    w.write_all(chunk).map_err(|e| format!("Write error: {e}"))?;
                    writeln!(w).map_err(|e| format!("Write error: {e}"))?;
                }
                count += 1;
            }
        }
    }

    w.flush().map_err(|e| format!("Flush error: {e}"))?;
    Ok(count)
}
