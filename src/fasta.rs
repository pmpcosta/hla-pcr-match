//! Fast FASTA parser (plain text and .zip) and allele query filtering for IMGT/HLA and generic MHC FASTAs.

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct AlleleRecord {
    /// Accession ID if present (e.g. "HLA:HLA00001")
    pub accession: String,
    /// Allele name (e.g. "A*01:01:01:01", "MICA*001:01:01", "Mamu-A1*001:01")
    pub name: String,
    /// Uppercase allele name cached for fast case-insensitive filtering
    pub name_upper: String,
    /// Locus extracted before '*' (e.g. "A", "DRB1", "MICA", "Mamu-A1")
    pub locus: String,
    /// Uppercase locus cached for fast case-insensitive filtering
    pub locus_upper: String,
    /// Uppercase ASCII nucleotide sequence
    pub sequence: Vec<u8>,
}

#[derive(Debug, Clone, Default)]
pub struct AlleleDatabase {
    pub source_label: String,
    pub records: Vec<AlleleRecord>,
    pub loci: Vec<String>,
}

impl AlleleDatabase {
    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Loads alleles from a plain FASTA file or a ZIP archive containing a FASTA file.
    pub fn load_from_path(path: &Path, label: String) -> Result<Self, String> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();

        if ext == "zip" {
            let file =
                File::open(path).map_err(|e| format!("Failed to open {}: {e}", path.display()))?;
            let mut archive = zip::ZipArchive::new(file)
                .map_err(|e| format!("Failed to read zip archive {}: {e}", path.display()))?;
            for i in 0..archive.len() {
                let entry = archive
                    .by_index(i)
                    .map_err(|e| format!("Failed to read zip entry {i}: {e}"))?;
                let name = entry.name().to_ascii_lowercase();
                if name.ends_with(".fasta") || name.ends_with(".fa") || name.ends_with(".fna") {
                    return Self::parse_reader(BufReader::new(entry), label);
                }
            }
            Err(format!(
                "No .fasta/.fa/.fna file found inside {}",
                path.display()
            ))
        } else {
            let file =
                File::open(path).map_err(|e| format!("Failed to open {}: {e}", path.display()))?;
            Self::parse_reader(BufReader::with_capacity(256 * 1024, file), label)
        }
    }

    pub fn parse_reader<R: Read>(reader: BufReader<R>, label: String) -> Result<Self, String> {
        let mut records = Vec::with_capacity(4096);
        let mut loci_set = BTreeSet::new();

        let mut current_header: Option<String> = None;
        let mut current_seq: Vec<u8> = Vec::with_capacity(4096);

        let mut line = String::new();
        let mut reader = reader;
        loop {
            line.clear();
            let bytes = reader
                .read_line(&mut line)
                .map_err(|e| format!("Error reading FASTA: {e}"))?;
            if bytes == 0 {
                break;
            }
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if let Some(header_body) = trimmed.strip_prefix('>') {
                if let Some(prev_header) = current_header.take() {
                    if !current_seq.is_empty() {
                        let rec = build_record(&prev_header, std::mem::take(&mut current_seq));
                        loci_set.insert(rec.locus.clone());
                        records.push(rec);
                    }
                }
                current_header = Some(header_body.to_string());
                current_seq.clear();
            } else if current_header.is_some() {
                for &b in trimmed.as_bytes() {
                    if !b.is_ascii_whitespace() {
                        current_seq.push(b.to_ascii_uppercase());
                    }
                }
            }
        }

        if let Some(prev_header) = current_header {
            if !current_seq.is_empty() {
                let rec = build_record(&prev_header, current_seq);
                loci_set.insert(rec.locus.clone());
                records.push(rec);
            }
        }

        if records.is_empty() {
            return Err("No valid FASTA sequences found in file".to_string());
        }

        Ok(Self {
            source_label: label,
            records,
            loci: loci_set.into_iter().collect(),
        })
    }
}

fn build_record(header: &str, sequence: Vec<u8>) -> AlleleRecord {
    let mut parts = header.split_whitespace();
    let first = parts.next().unwrap_or("");
    let second = parts.next();

    let (accession, name) = if (first.starts_with("HLA:") || first.starts_with("IPD:"))
        && second.is_some()
    {
        (first.to_string(), second.unwrap().to_string())
    } else {
        (String::new(), first.to_string())
    };

    let locus = match name.split_once('*') {
        Some((loc, _)) => loc.to_string(),
        None => name.clone(),
    };

    let name_upper = name.to_ascii_uppercase();
    let locus_upper = locus.to_ascii_uppercase();

    AlleleRecord {
        accession,
        name,
        name_upper,
        locus,
        locus_upper,
        sequence,
    }
}

#[derive(Debug, Clone)]
pub struct QueryToken {
    pub raw: String,
    normalized: String,
    stripped_hla: Option<String>,
}

impl QueryToken {
    fn from_raw(raw: &str) -> Self {
        let upper = raw.trim().to_ascii_uppercase();
        let stripped_hla = upper.strip_prefix("HLA-").map(|s| s.to_string());
        Self {
            raw: raw.trim().to_string(),
            normalized: upper,
            stripped_hla,
        }
    }

    pub fn matches(&self, rec: &AlleleRecord) -> bool {
        if token_matches_record(&self.normalized, rec) {
            return true;
        }
        if let Some(ref stripped) = self.stripped_hla {
            if token_matches_record(stripped, rec) {
                return true;
            }
        }
        false
    }
}

fn token_matches_record(token: &str, rec: &AlleleRecord) -> bool {
    if token.is_empty() {
        return false;
    }

    // e.g. "A*" or "DRB1*" -> match entire locus
    if let Some(loc_only) = token.strip_suffix('*') {
        if !loc_only.contains('*') {
            return rec.locus_upper == loc_only;
        }
        return rec.name_upper.starts_with(loc_only);
    }

    if token.contains('*') {
        // Specific allele or allele-group prefix (e.g. "A*02", "A*02:01", "A*02:01:01:01")
        if rec.name_upper == token {
            return true;
        }
        if let Some(rest) = rec.name_upper.strip_prefix(token) {
            if let Some(next_ch) = rest.chars().next() {
                // Respect field boundaries (':') or expression suffixes (N, L, S, C, A, Q)
                return next_ch == ':' || next_ch.is_ascii_alphabetic();
            }
        }
        false
    } else {
        // Locus-only query (e.g. "A", "DRB1", "MICA", "Mamu-A1") or accession match
        rec.locus_upper == token || rec.name_upper == token
    }
}

/// Splits the query textarea by newlines, commas, or semicolons.
pub fn parse_query_tokens(query: &str) -> Vec<QueryToken> {
    query
        .split(|c| c == ',' || c == '\n' || c == '\r' || c == ';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(QueryToken::from_raw)
        .collect()
}

/// Selects matching allele indices and returns any query tokens that matched zero alleles.
pub fn filter_alleles(db: &AlleleDatabase, query: &str) -> (Vec<usize>, Vec<String>) {
    let tokens = parse_query_tokens(query);
    if tokens.is_empty() {
        return ((0..db.records.len()).collect(), Vec::new());
    }

    let mut token_hit_counts = vec![0usize; tokens.len()];
    let mut selected = Vec::new();

    for (idx, rec) in db.records.iter().enumerate() {
        let mut matched = false;
        for (ti, tok) in tokens.iter().enumerate() {
            if tok.matches(rec) {
                token_hit_counts[ti] += 1;
                matched = true;
            }
        }
        if matched {
            selected.push(idx);
        }
    }

    let unmatched_queries = tokens
        .into_iter()
        .zip(token_hit_counts)
        .filter_map(|(tok, cnt)| if cnt == 0 { Some(tok.raw) } else { None })
        .collect();

    (selected, unmatched_queries)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SAMPLE_FASTA: &str = "\
>HLA:HLA00001 A*01:01:01:01 20 bp
ATGGCCGTCATGGCGCCCCG
>HLA:HLA00002 A*01:02 20 bp
ATGGCCGTCATGGCGCCCCA
>HLA:HLA00003 A*010:01 20 bp
ATGGCCGTCATGGCGCCCCT
>HLA:HLA00100 B*07:02:01:01 20 bp
ATGCGGGTCACGGCGCCCCG
>HLA:HLA01013 MICA*001:01:01 20 bp
ATGGGGCTGGGCCCGGTCTT
";

    #[test]
    fn test_parse_and_filter_hla_queries() {
        let db = AlleleDatabase::parse_reader(
            BufReader::new(Cursor::new(SAMPLE_FASTA.as_bytes())),
            "test".to_string(),
        )
        .unwrap();
        assert_eq!(db.len(), 5);
        assert_eq!(db.loci, vec!["A", "B", "MICA"]);

        // HLA-A matches locus A
        let (ids, unmatched) = filter_alleles(&db, "HLA-A");
        assert_eq!(ids, vec![0, 1, 2]);
        assert!(unmatched.is_empty());

        // HLA-A*01 respects field boundary and does NOT match A*010:01
        let (ids, unmatched) = filter_alleles(&db, "HLA-A*01");
        assert_eq!(ids, vec![0, 1]);
        assert!(unmatched.is_empty());

        // Multi-token + case insensitive + unmatched reporting
        let (ids, unmatched) = filter_alleles(&db, "hla-b*07:02, mica\nHLA-Z99");
        assert_eq!(ids, vec![3, 4]);
        assert_eq!(unmatched, vec!["HLA-Z99"]);

        // Empty query selects all
        let (ids, unmatched) = filter_alleles(&db, "   ");
        assert_eq!(ids.len(), 5);
        assert!(unmatched.is_empty());
    }
}
