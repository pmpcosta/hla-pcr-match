//! Parallel in-silico PCR primer matching engine modeled after IPD-MHC Primer Match.

use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

use rayon::prelude::*;

use crate::fasta::{AlleleDatabase, filter_alleles};
use crate::iupac::{bases_match, clean_and_validate_primer, encode_primer, reverse_complement};

/// Safety limit on amplicon combinations recorded per allele to prevent OOM on
/// highly degenerate or repetitive primers.
pub const MAX_AMPLICONS_PER_ALLELE: usize = 64;

#[derive(Debug, Clone)]
pub struct MatchParams {
    pub query: String,
    pub primer_fwd: String,
    pub primer_rvs: String,
    /// Number of bases at the 3' end of each primer that must match identically (0 mismatches).
    pub end_identity: usize,
    /// Maximum total mismatches allowed across the full primer sequence.
    pub max_mismatch: usize,
    /// If true, reverse-complement the reverse primer before scanning the sense strand.
    pub flip_reverse: bool,
    /// Minimum amplicon size in bp (-1 = no lower bound).
    pub product_size_min: i64,
    /// Maximum amplicon size in bp (-1 = no upper bound).
    pub product_size_max: i64,
}

impl Default for MatchParams {
    fn default() -> Self {
        Self {
            query: String::new(),
            primer_fwd: String::new(),
            primer_rvs: String::new(),
            end_identity: 3,
            max_mismatch: 3,
            flip_reverse: true,
            product_size_min: -1,
            product_size_max: -1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct AmpliconHit {
    pub allele_index: usize,
    pub allele: String,
    pub accession: String,
    pub locus: String,
    /// 0-based inclusive start coordinate on template (forward primer 5' start).
    pub start: usize,
    /// 0-based exclusive end coordinate on template (reverse primer 5' end on sense strand).
    pub end: usize,
    /// Amplicon length (`end - start`).
    pub length: usize,
    pub fwd_mismatches: usize,
    pub rev_mismatches: usize,
    /// Template sequence at forward primer binding site (5'->3' sense strand).
    pub fwd_template: String,
    /// Template sequence at reverse primer binding site (5'->3' sense strand).
    pub rev_template: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MissReason {
    NeitherPrimerBound,
    ForwardNotBound,
    ReverseNotBound,
    InvalidOrientation,
    ProductSizeOutOfRange { min_found: usize, max_found: usize },
}

impl MissReason {
    pub fn label(&self) -> String {
        match self {
            Self::NeitherPrimerBound => "Neither primer bound".to_string(),
            Self::ForwardNotBound => "Forward primer did not bind".to_string(),
            Self::ReverseNotBound => "Reverse primer did not bind".to_string(),
            Self::InvalidOrientation => "Reverse site upstream of forward site".to_string(),
            Self::ProductSizeOutOfRange {
                min_found,
                max_found,
            } => {
                if min_found == max_found {
                    format!("Product size out of range ({min_found} bp)")
                } else {
                    format!("Product size out of range ({min_found}–{max_found} bp)")
                }
            }
        }
    }

    /// Zero-allocation sort key matching the alphabetical order of `label()`.
    pub fn sort_key(&self) -> (u8, usize, usize) {
        match self {
            Self::ForwardNotBound => (0, 0, 0),
            Self::NeitherPrimerBound => (1, 0, 0),
            Self::ProductSizeOutOfRange {
                min_found,
                max_found,
            } => (2, *min_found, *max_found),
            Self::ReverseNotBound => (3, 0, 0),
            Self::InvalidOrientation => (4, 0, 0),
        }
    }
}

#[derive(Debug, Clone)]
pub struct MissedAllele {
    pub allele: String,
    pub accession: String,
    pub locus: String,
    pub reason: MissReason,
    pub fwd_hits: usize,
    pub rev_hits: usize,
}

#[derive(Debug, Clone)]
pub struct LocusCoverage {
    pub locus: String,
    pub total_alleles: usize,
    pub matched_alleles: usize,
    pub missed_alleles: usize,
    pub coverage_pct: f64,
}

#[derive(Debug, Clone)]
pub struct MatchResults {
    /// Database snapshot these results were computed against (keeps `allele_index` valid).
    pub db: Arc<AlleleDatabase>,
    pub hits: Vec<AmpliconHit>,
    pub missed: Vec<MissedAllele>,
    pub locus_coverage: Vec<LocusCoverage>,
    pub total_searched: usize,
    pub matched_unique_alleles: usize,
    pub elapsed_ms: u128,
    pub unmatched_queries: Vec<String>,
    /// True if at least one allele hit `MAX_AMPLICONS_PER_ALLELE` and was truncated.
    pub truncated_hits: bool,
    /// Cleaned forward primer (5'->3' sense strand).
    pub fwd_primer_used: String,
    /// Effective reverse primer motif scanned on the sense strand (RC if flip_reverse=true).
    pub rev_motif_used: String,
    /// End identity used for alignment highlighting.
    pub end_identity_used: usize,
}

#[derive(Debug, Clone, Copy)]
struct SiteHit {
    pos: usize,
    mismatches: usize,
}

/// Scans `seq` for forward primer matches where the LAST `end_id` bases must match exactly.
fn find_forward_hits(
    seq: &[u8],
    primer: &[u8],
    end_identity: usize,
    max_mismatch: usize,
) -> Vec<SiteHit> {
    let plen = primer.len();
    if seq.len() < plen || plen == 0 {
        return Vec::new();
    }
    let end_id = end_identity.min(plen);
    let prefix_len = plen - end_id;
    let mut hits = Vec::new();

    'outer: for pos in 0..=(seq.len() - plen) {
        let window = &seq[pos..pos + plen];

        // Fast fail on 3'-end exact identity region (last `end_id` bases)
        for i in prefix_len..plen {
            if !bases_match(window[i], primer[i]) {
                continue 'outer;
            }
        }

        // Count mismatches over remaining 5' portion with early exit
        let mut mm = 0usize;
        for i in 0..prefix_len {
            if !bases_match(window[i], primer[i]) {
                mm += 1;
                if mm > max_mismatch {
                    continue 'outer;
                }
            }
        }

        hits.push(SiteHit {
            pos,
            mismatches: mm,
        });
    }
    hits
}

/// Scans `seq` for reverse primer matches (already oriented 5'->3' on sense strand),
/// where the FIRST `end_id` bases correspond to the reverse primer's 3' end and must match exactly.
fn find_reverse_hits(
    seq: &[u8],
    rev_sense_motif: &[u8],
    end_identity: usize,
    max_mismatch: usize,
) -> Vec<SiteHit> {
    let plen = rev_sense_motif.len();
    if seq.len() < plen || plen == 0 {
        return Vec::new();
    }
    let end_id = end_identity.min(plen);
    let mut hits = Vec::new();

    'outer: for pos in 0..=(seq.len() - plen) {
        let window = &seq[pos..pos + plen];

        // Fast fail on reverse primer's 3'-end exact identity region (first `end_id` bases on sense strand)
        for i in 0..end_id {
            if !bases_match(window[i], rev_sense_motif[i]) {
                continue 'outer;
            }
        }

        // Count mismatches over remaining bases with early exit
        let mut mm = 0usize;
        for i in end_id..plen {
            if !bases_match(window[i], rev_sense_motif[i]) {
                mm += 1;
                if mm > max_mismatch {
                    continue 'outer;
                }
            }
        }

        hits.push(SiteHit {
            pos,
            mismatches: mm,
        });
    }
    hits
}

enum AlleleOutcome {
    Matched {
        hits: Vec<AmpliconHit>,
        truncated: bool,
    },
    Missed(MissedAllele),
}

pub fn run_match(db: Arc<AlleleDatabase>, params: &MatchParams) -> Result<MatchResults, String> {
    let started = Instant::now();

    if db.is_empty() {
        return Err(
            "No allele database loaded. Download IMGT/HLA or load a FASTA file first.".into(),
        );
    }

    let fwd_clean = clean_and_validate_primer(&params.primer_fwd, "Forward primer")?;
    let rvs_clean = clean_and_validate_primer(&params.primer_rvs, "Reverse primer")?;

    if params.product_size_min >= 0
        && params.product_size_max >= 0
        && params.product_size_min > params.product_size_max
    {
        return Err("Product Size Min cannot be greater than Product Size Max.".into());
    }

    let rev_sense_str = if params.flip_reverse {
        reverse_complement(&rvs_clean)
    } else {
        rvs_clean.clone()
    };

    let fwd_mask = encode_primer(&fwd_clean);
    let rev_mask = encode_primer(&rev_sense_str);
    let fwd_len = fwd_mask.len();
    let rev_len = rev_mask.len();

    let (selected_indices, unmatched_queries) = filter_alleles(&db, &params.query);
    let total_searched = selected_indices.len();

    if total_searched == 0 {
        return Err(format!(
            "No alleles matched query '{}'. Try e.g. HLA-A, DRB1, A*02:01, or leave empty to search all loci.",
            params.query.trim()
        ));
    }

    let min_size = if params.product_size_min >= 0 {
        Some(params.product_size_min as usize)
    } else {
        None
    };
    let max_size = if params.product_size_max >= 0 {
        Some(params.product_size_max as usize)
    } else {
        None
    };

    let records = &db.records;
    let outcomes: Vec<AlleleOutcome> = selected_indices
        .par_iter()
        .map(|&allele_idx| {
            let rec = &records[allele_idx];
            let seq = &rec.sequence;

            let f_hits =
                find_forward_hits(seq, &fwd_mask, params.end_identity, params.max_mismatch);
            let r_hits =
                find_reverse_hits(seq, &rev_mask, params.end_identity, params.max_mismatch);

            if f_hits.is_empty() || r_hits.is_empty() {
                let reason = match (f_hits.is_empty(), r_hits.is_empty()) {
                    (true, true) => MissReason::NeitherPrimerBound,
                    (true, false) => MissReason::ForwardNotBound,
                    (false, true) => MissReason::ReverseNotBound,
                    (false, false) => unreachable!(),
                };
                return AlleleOutcome::Missed(MissedAllele {
                    allele: rec.name.clone(),
                    accession: rec.accession.clone(),
                    locus: rec.locus.clone(),
                    reason,
                    fwd_hits: f_hits.len(),
                    rev_hits: r_hits.len(),
                });
            }

            let mut allele_amplicons = Vec::new();
            let mut valid_orientation_min: Option<usize> = None;
            let mut valid_orientation_max: Option<usize> = None;
            let mut truncated = false;

            for fh in &f_hits {
                // Amplicon must start strictly before reverse site and extend strictly
                // past the 3' end of the forward primer (rejects reversed and nested sites).
                let min_r_pos = (fh.pos + 1).max((fh.pos + fwd_len + 1).saturating_sub(rev_len));
                let r_start = r_hits.partition_point(|rh| rh.pos < min_r_pos);
                if r_start >= r_hits.len() {
                    continue;
                }

                let first_len = r_hits[r_start].pos + rev_len - fh.pos;
                let last_len = r_hits[r_hits.len() - 1].pos + rev_len - fh.pos;
                valid_orientation_min =
                    Some(valid_orientation_min.map_or(first_len, |m| m.min(first_len)));
                valid_orientation_max =
                    Some(valid_orientation_max.map_or(last_len, |m| m.max(last_len)));

                let r_min_idx = match min_size {
                    Some(min_s) => {
                        let target_pos = fh.pos.saturating_add(min_s).saturating_sub(rev_len);
                        r_start + r_hits[r_start..].partition_point(|rh| rh.pos < target_pos)
                    }
                    None => r_start,
                };

                for rh in &r_hits[r_min_idx..] {
                    let end = rh.pos + rev_len;
                    let length = end - fh.pos;
                    if max_size.is_some_and(|max_s| length > max_s) {
                        break;
                    }
                    if allele_amplicons.len() >= MAX_AMPLICONS_PER_ALLELE {
                        truncated = true;
                        break;
                    }

                    let fwd_template =
                        String::from_utf8_lossy(&seq[fh.pos..fh.pos + fwd_len]).into_owned();
                    let rev_template =
                        String::from_utf8_lossy(&seq[rh.pos..rh.pos + rev_len]).into_owned();
                    allele_amplicons.push(AmpliconHit {
                        allele_index: allele_idx,
                        allele: rec.name.clone(),
                        accession: rec.accession.clone(),
                        locus: rec.locus.clone(),
                        start: fh.pos,
                        end,
                        length,
                        fwd_mismatches: fh.mismatches,
                        rev_mismatches: rh.mismatches,
                        fwd_template,
                        rev_template,
                    });
                }

                if truncated {
                    break;
                }
            }

            if allele_amplicons.is_empty() {
                let reason = match (valid_orientation_min, valid_orientation_max) {
                    (Some(min_found), Some(max_found)) => MissReason::ProductSizeOutOfRange {
                        min_found,
                        max_found,
                    },
                    _ => MissReason::InvalidOrientation,
                };
                AlleleOutcome::Missed(MissedAllele {
                    allele: rec.name.clone(),
                    accession: rec.accession.clone(),
                    locus: rec.locus.clone(),
                    reason,
                    fwd_hits: f_hits.len(),
                    rev_hits: r_hits.len(),
                })
            } else {
                AlleleOutcome::Matched {
                    hits: allele_amplicons,
                    truncated,
                }
            }
        })
        .collect();

    let mut hits = Vec::new();
    let mut missed = Vec::new();
    let mut matched_unique_alleles = 0usize;
    let mut truncated_hits = false;
    let mut locus_totals: BTreeMap<String, (usize, usize)> = BTreeMap::new(); // (matched, missed)

    for outcome in outcomes {
        match outcome {
            AlleleOutcome::Matched {
                hits: mut allele_hits,
                truncated,
            } => {
                matched_unique_alleles += 1;
                truncated_hits |= truncated;
                if let Some(first) = allele_hits.first() {
                    locus_totals.entry(first.locus.clone()).or_default().0 += 1;
                }
                hits.append(&mut allele_hits);
            }
            AlleleOutcome::Missed(miss) => {
                locus_totals.entry(miss.locus.clone()).or_default().1 += 1;
                missed.push(miss);
            }
        }
    }

    let locus_coverage: Vec<LocusCoverage> = locus_totals
        .into_iter()
        .map(|(locus, (matched, missed_cnt))| {
            let total = matched + missed_cnt;
            let pct = if total > 0 {
                (matched as f64 / total as f64) * 100.0
            } else {
                0.0
            };
            LocusCoverage {
                locus,
                total_alleles: total,
                matched_alleles: matched,
                missed_alleles: missed_cnt,
                coverage_pct: pct,
            }
        })
        .collect();

    Ok(MatchResults {
        db,
        hits,
        missed,
        locus_coverage,
        total_searched,
        matched_unique_alleles,
        elapsed_ms: started.elapsed().as_millis(),
        unmatched_queries,
        truncated_hits,
        fwd_primer_used: fwd_clean,
        rev_motif_used: rev_sense_str,
        end_identity_used: params.end_identity,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    // Real HLA-A*01:01:01:01 prefix (first 258 bp from ANHIG/IMGTHLA A_nuc.fasta)
    const TEST_FASTA: &str = "\
>HLA:HLA00001 A*01:01:01:01 258 bp
ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAG
>HLA:HLA00002 A*01:02 258 bp
ATGGCCGTCATGGCGCCCCTAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAG
>HLA:HLA00003 A*02:01:01:01 258 bp
ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACGGG
";

    fn parse_db(fasta: &str) -> Arc<AlleleDatabase> {
        Arc::new(
            AlleleDatabase::parse_reader(
                BufReader::new(Cursor::new(fasta.as_bytes())),
                "test".to_string(),
            )
            .unwrap(),
        )
    }

    #[test]
    fn test_exact_and_end_identity_matching() {
        let db = parse_db(TEST_FASTA);

        // Forward primer matches 0..20 of A*01:01:01:01 exactly; A*01:02 has a mismatch at base 19 (3' end!)
        // Reverse primer (5'->3') is RC of last 18 bp ("CCGGAGTATTGGGACCAG" -> "CTGGTCCCAATACTCCGG")
        let params = MatchParams {
            query: "HLA-A".to_string(),
            primer_fwd: "ATGGCCGTCATGGCGCCCCG".to_string(),
            primer_rvs: "CTGGTCCCAATACTCCGG".to_string(), // RC = CCGGAGTATTGGGACCAG (pos 240..258)
            end_identity: 3,
            max_mismatch: 3,
            flip_reverse: true,
            product_size_min: -1,
            product_size_max: -1,
        };

        let res = run_match(db, &params).unwrap();
        // A*01:01:01:01 matches exact (0..258, len 258).
        // A*01:02 fails forward 3'-end identity (base 19 is T vs G).
        // A*02:01:01:01 has "GGG" at 255..258, which is the 5' end of the reverse primer (2 mismatches, allowed).
        assert_eq!(res.matched_unique_alleles, 2);
        assert_eq!(res.hits[0].allele, "A*01:01:01:01");
        assert_eq!(res.hits[0].start, 0);
        assert_eq!(res.hits[0].end, 258);
        assert_eq!(res.hits[0].length, 258);
        assert_eq!(res.hits[0].fwd_mismatches, 0);
        assert_eq!(res.hits[0].rev_mismatches, 0);

        assert_eq!(res.hits[1].allele, "A*02:01:01:01");
        assert_eq!(res.hits[1].rev_mismatches, 2);

        assert_eq!(res.missed.len(), 1);
        assert_eq!(res.missed[0].allele, "A*01:02");
        assert_eq!(res.missed[0].reason, MissReason::ForwardNotBound);
    }

    #[test]
    fn test_reverse_3prime_end_identity_rejected() {
        let db = parse_db(TEST_FASTA);

        // Reverse primer 5'-CTGGTCCCAATACTCCGT-3' -> 3' base is T -> RC 5' base is A (template has C -> violates rev 3' end identity!)
        let params = MatchParams {
            query: "A*01:01:01:01".to_string(),
            primer_fwd: "ATGGCCGTCATGGCGCCCCG".to_string(),
            primer_rvs: "CTGGTCCCAATACTCCGT".to_string(),
            end_identity: 3,
            max_mismatch: 3,
            flip_reverse: true,
            product_size_min: -1,
            product_size_max: -1,
        };

        let res = run_match(Arc::clone(&db), &params).unwrap();
        assert_eq!(res.hits.len(), 0);
        assert_eq!(res.missed.len(), 1);
        assert_eq!(res.missed[0].reason, MissReason::ReverseNotBound);

        // With end_identity = 0, that 1 mismatch IS allowed!
        let mut params_no_ei = params;
        params_no_ei.end_identity = 0;
        let res2 = run_match(db, &params_no_ei).unwrap();
        assert_eq!(res2.hits.len(), 1);
        assert_eq!(res2.hits[0].rev_mismatches, 1);
    }

    #[test]
    fn test_product_size_filter() {
        let db = parse_db(TEST_FASTA);

        let params = MatchParams {
            query: "A*01:01:01:01".to_string(),
            primer_fwd: "ATGGCCGTCATGGCGCCCCG".to_string(),
            primer_rvs: "CTGGTCCCAATACTCCGG".to_string(),
            end_identity: 3,
            max_mismatch: 3,
            flip_reverse: true,
            product_size_min: 300,
            product_size_max: 500,
        };

        let res = run_match(db, &params).unwrap();
        assert_eq!(res.hits.len(), 0);
        assert_eq!(
            res.missed[0].reason,
            MissReason::ProductSizeOutOfRange {
                min_found: 258,
                max_found: 258
            }
        );
    }

    #[test]
    fn test_flip_reverse_false_and_degenerate_primer() {
        let db = parse_db(TEST_FASTA);

        // Provide reverse primer already oriented on the sense strand, with degenerate Y matching C/T
        let params = MatchParams {
            query: "A*01:01:01:01".to_string(),
            primer_fwd: "ATGGCCGTCATGGCGCCCYG".to_string(),
            primer_rvs: "CCGGAGTATTGGGACCAG".to_string(),
            end_identity: 3,
            max_mismatch: 0,
            flip_reverse: false,
            product_size_min: 200,
            product_size_max: 300,
        };

        let res = run_match(db, &params).unwrap();
        assert_eq!(res.hits.len(), 1);
        assert_eq!(res.hits[0].length, 258);
        assert_eq!(res.hits[0].fwd_mismatches, 0);
        assert_eq!(res.hits[0].rev_mismatches, 0);
    }

    #[test]
    fn test_invalid_orientation_and_nested_reverse_rejected() {
        // Reverse site is nested inside forward primer (fwd len 12 at 0..12, rev len 6 at 2..8 and 6..12)
        let fasta = ">HLA:HLA99999 X*01:01\nACGTACGTACGTGGGGGG\n";
        let db = parse_db(fasta);

        let params = MatchParams {
            query: String::new(),
            primer_fwd: "ACGTACGTACGT".to_string(),
            primer_rvs: "GTACGT".to_string(), // sense motif directly (flip_reverse=false), matches 2..8 and 6..12
            end_identity: 2,
            max_mismatch: 0,
            flip_reverse: false,
            product_size_min: -1,
            product_size_max: -1,
        };

        let res = run_match(db, &params).unwrap();
        assert_eq!(res.hits.len(), 0);
        assert_eq!(res.missed.len(), 1);
        assert_eq!(res.missed[0].reason, MissReason::InvalidOrientation);

        // Reverse site strictly upstream of forward site
        let fasta2 = ">HLA:HLA99998 X*01:02\nCCCCCCTTTTTTGGGGGGGG\n";
        let db2 = parse_db(fasta2);
        let params2 = MatchParams {
            query: String::new(),
            primer_fwd: "GGGGGGGG".to_string(),
            primer_rvs: "CCCCCC".to_string(),
            end_identity: 2,
            max_mismatch: 0,
            flip_reverse: false,
            product_size_min: -1,
            product_size_max: -1,
        };
        let res2 = run_match(db2, &params2).unwrap();
        assert_eq!(res2.hits.len(), 0);
        assert_eq!(res2.missed[0].reason, MissReason::InvalidOrientation);
    }

    #[test]
    fn test_poly_n_gap_does_not_match() {
        let fasta = ">HLA:HLA88888 X*02:01\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";
        let db = parse_db(fasta);

        let params = MatchParams {
            query: String::new(),
            primer_fwd: "ACGTACGT".to_string(),
            primer_rvs: "NNNNNNNN".to_string(),
            end_identity: 0,
            max_mismatch: 2,
            flip_reverse: false,
            product_size_min: -1,
            product_size_max: -1,
        };

        let res = run_match(db, &params).unwrap();
        assert_eq!(res.hits.len(), 0);
        assert_eq!(res.missed.len(), 1);
        assert_eq!(res.missed[0].reason, MissReason::NeitherPrimerBound);
    }
}
