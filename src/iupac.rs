//! Bitmask-based IUPAC DNA nucleotide encoding, validation, and reverse complement.

/// 4-bit nucleotide bitmask table indexed by ASCII byte.
/// A = 0b0001, C = 0b0010, G = 0b0100, T/U = 0b1000.
/// Invalid/gap characters map to 0 (never matches any base).
pub const BASE_MASK: [u8; 256] = {
    let mut table = [0u8; 256];
    const A: u8 = 0b0001;
    const C: u8 = 0b0010;
    const G: u8 = 0b0100;
    const T: u8 = 0b1000;

    let pairs: &[(u8, u8)] = &[
        (b'A', A),
        (b'C', C),
        (b'G', G),
        (b'T', T),
        (b'U', T),
        (b'R', A | G),
        (b'Y', C | T),
        (b'S', G | C),
        (b'W', A | T),
        (b'K', G | T),
        (b'M', A | C),
        (b'B', C | G | T),
        (b'D', A | G | T),
        (b'H', A | C | T),
        (b'V', A | C | G),
        (b'N', A | C | G | T),
    ];

    let mut i = 0;
    while i < pairs.len() {
        let (ch, mask) = pairs[i];
        table[ch as usize] = mask;
        table[(ch + 32) as usize] = mask; // lowercase
        i += 1;
    }
    table
};

/// Returns true if the template base and primer base match under IUPAC degeneracy.
#[inline(always)]
pub fn bases_match(seq_byte: u8, primer_mask: u8) -> bool {
    (BASE_MASK[seq_byte as usize] & primer_mask) != 0
}

/// Cleans whitespace, normalizes to uppercase (U -> T), and validates all bases are IUPAC DNA.
pub fn clean_and_validate_primer(raw: &str, label: &str) -> Result<String, String> {
    let mut cleaned = String::with_capacity(raw.len());
    for ch in raw.chars() {
        if ch.is_whitespace() {
            continue;
        }
        let up = ch.to_ascii_uppercase();
        let norm = if up == 'U' { 'T' } else { up };
        if !norm.is_ascii() || BASE_MASK[norm as usize] == 0 {
            return Err(format!(
                "{label} contains invalid IUPAC nucleotide character '{ch}'"
            ));
        }
        cleaned.push(norm);
    }
    if cleaned.is_empty() {
        return Err(format!("{label} cannot be empty"));
    }
    Ok(cleaned)
}

/// Encodes an already-validated primer string into 4-bit IUPAC masks.
pub fn encode_primer(primer: &str) -> Vec<u8> {
    primer.bytes().map(|b| BASE_MASK[b as usize]).collect()
}

/// Computes the IUPAC reverse complement of a validated uppercase DNA string.
pub fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'T' | b'U' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            b'R' => 'Y',
            b'Y' => 'R',
            b'K' => 'M',
            b'M' => 'K',
            b'S' => 'S',
            b'W' => 'W',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            b'N' => 'N',
            other => other as char,
        })
        .collect()
}

/// Calculates total sequence degeneracy (product of base possibilities per position).
pub fn primer_degeneracy(primer: &str) -> u64 {
    let mut total: u64 = 1;
    for b in primer.bytes() {
        let m = BASE_MASK[b as usize];
        if m != 0 {
            total = total.saturating_mul(m.count_ones() as u64);
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_iupac() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("ATGGCRYSWKMBDHVN"), "NBDHVKMWSRYGCCAT");
        assert_eq!(reverse_complement("TGCGAGTCCCAATACTCTGG"), "CCAGAGTATTGGGACTCGCA");
    }

    #[test]
    fn test_clean_and_validate() {
        assert_eq!(
            clean_and_validate_primer(" atg gcu r ", "Fwd").unwrap(),
            "ATGGCTR"
        );
        assert!(clean_and_validate_primer("ATGX", "Fwd").is_err());
        assert!(clean_and_validate_primer("   ", "Fwd").is_err());
    }

    #[test]
    fn test_degeneracy_and_matching() {
        assert_eq!(primer_degeneracy("ACGT"), 1);
        assert_eq!(primer_degeneracy("ACRYN"), 16);
        assert!(bases_match(b'A', BASE_MASK[b'R' as usize]));
        assert!(bases_match(b'G', BASE_MASK[b'R' as usize]));
        assert!(!bases_match(b'C', BASE_MASK[b'R' as usize]));
        assert!(!bases_match(b'-', BASE_MASK[b'N' as usize]));
    }
}
