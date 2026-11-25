//! CSV, TSV, coverage summary, and amplicon FASTA export utilities.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::matcher::MatchResults;

fn escape_field(value: &str, sep: char) -> String {
    if value.contains(sep) || value.contains('"') || value.contains('\n') || value.contains('\r') {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

fn write_hits_delimited<W: Write>(
    w: &mut W,
    results: &MatchResults,
    sep: char,
) -> Result<usize, String> {
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

pub fn export_hits_delimited(
    path: &Path,
    results: &MatchResults,
    sep: char,
) -> Result<usize, String> {
    let file =
        File::create(path).map_err(|e| format!("Failed to create {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);
    write_hits_delimited(&mut w, results, sep)
}

fn write_missed_csv<W: Write>(w: &mut W, results: &MatchResults) -> Result<usize, String> {
    writeln!(
        w,
        "Allele,Accession,Locus,Reason,FwdBindingSites,RevBindingSites"
    )
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

pub fn export_missed_csv(path: &Path, results: &MatchResults) -> Result<usize, String> {
    let file =
        File::create(path).map_err(|e| format!("Failed to create {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);
    write_missed_csv(&mut w, results)
}

fn write_coverage_csv<W: Write>(w: &mut W, results: &MatchResults) -> Result<usize, String> {
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

pub fn export_coverage_csv(path: &Path, results: &MatchResults) -> Result<usize, String> {
    let file =
        File::create(path).map_err(|e| format!("Failed to create {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);
    write_coverage_csv(&mut w, results)
}

fn write_amplicons_fasta<W: Write>(w: &mut W, results: &MatchResults) -> Result<usize, String> {
    let mut count = 0usize;
    for h in &results.hits {
        if let Some(rec) = results.db.records.get(h.allele_index)
            && h.end <= rec.sequence.len()
            && h.start < h.end
        {
            let amplicon = &rec.sequence[h.start..h.end];
            if h.accession.is_empty() {
                writeln!(
                    w,
                    ">{} start={} end={} len={} fwd_mm={} rev_mm={}",
                    h.allele, h.start, h.end, h.length, h.fwd_mismatches, h.rev_mismatches
                )
                .map_err(|e| format!("Write error: {e}"))?;
            } else {
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
            }

            for chunk in amplicon.chunks(80) {
                w.write_all(chunk)
                    .map_err(|e| format!("Write error: {e}"))?;
                writeln!(w).map_err(|e| format!("Write error: {e}"))?;
            }
            count += 1;
        }
    }

    w.flush().map_err(|e| format!("Flush error: {e}"))?;
    Ok(count)
}

pub fn export_amplicons_fasta(path: &Path, results: &MatchResults) -> Result<usize, String> {
    let file =
        File::create(path).map_err(|e| format!("Failed to create {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);
    write_amplicons_fasta(&mut w, results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};
    use std::sync::Arc;

    use crate::fasta::AlleleDatabase;
    use crate::matcher::{MatchParams, run_match};

    fn sample_results() -> MatchResults {
        let fasta = "\
>HLA:HLA00001 A*01:01 30 bp
ATGGCCGTCACCCCCCCCCCCTGGTCCCAA
>CustomSeq
ATGGCCGTCACCCCCCCCCCCTGGTCCCAA
>HLA:HLA00002 A*01:02 30 bp
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
";
        let db = Arc::new(
            AlleleDatabase::parse_reader(
                BufReader::new(Cursor::new(fasta.as_bytes())),
                "test".to_string(),
            )
            .unwrap(),
        );
        let params = MatchParams {
            query: String::new(),
            primer_fwd: "ATGGCCGTCA".to_string(),
            primer_rvs: "TTGGGACCAG".to_string(), // RC = CTGGTCCCAA
            end_identity: 3,
            max_mismatch: 0,
            flip_reverse: true,
            product_size_min: -1,
            product_size_max: -1,
        };
        run_match(db, &params).unwrap()
    }

    #[test]
    fn test_escape_field() {
        assert_eq!(escape_field("plain", ','), "plain");
        assert_eq!(escape_field("a,b", ','), "\"a,b\"");
        assert_eq!(escape_field("a\"b", ','), "\"a\"\"b\"");
        assert_eq!(escape_field("a,b", '\t'), "a,b");
        assert_eq!(escape_field("a\tb", '\t'), "\"a\tb\"");
    }

    #[test]
    fn test_export_formats() {
        let res = sample_results();

        let mut csv_buf = Vec::new();
        let hit_cnt = write_hits_delimited(&mut csv_buf, &res, ',').unwrap();
        assert_eq!(hit_cnt, 2);
        let csv = String::from_utf8(csv_buf).unwrap();
        assert!(csv.starts_with("Allele,Accession,Locus,Start,End,Length,"));
        assert!(csv.contains("A*01:01,HLA:HLA00001,A,0,30,30,0,0,0,ATGGCCGTCA,CTGGTCCCAA"));

        let mut miss_buf = Vec::new();
        let miss_cnt = write_missed_csv(&mut miss_buf, &res).unwrap();
        assert_eq!(miss_cnt, 1);
        let miss = String::from_utf8(miss_buf).unwrap();
        assert!(miss.contains("A*01:02,HLA:HLA00002,A,Neither primer bound,0,0"));

        let mut cov_buf = Vec::new();
        let cov_cnt = write_coverage_csv(&mut cov_buf, &res).unwrap();
        assert_eq!(cov_cnt, 2);
        let cov = String::from_utf8(cov_buf).unwrap();
        assert!(cov.contains("A,1,1,2,50.00"));
        assert!(cov.contains("CustomSeq,1,0,1,100.00"));

        let mut fa_buf = Vec::new();
        let fa_cnt = write_amplicons_fasta(&mut fa_buf, &res).unwrap();
        assert_eq!(fa_cnt, 2);
        let fa = String::from_utf8(fa_buf).unwrap();
        assert!(fa.contains(">A*01:01 HLA:HLA00001 start=0 end=30 len=30 fwd_mm=0 rev_mm=0\n"));
        // Empty accession must not introduce a double space
        assert!(fa.contains(">CustomSeq start=0 end=30 len=30 fwd_mm=0 rev_mm=0\n"));
    }
}
