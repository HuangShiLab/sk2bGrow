//! Nucleotide primitives: IUPAC matching, reverse complement, 2-bit packing, GC.
//!
//! Mirrors the helpers in Syn2b's `bsyn::seq`; kept vendored so the workspace
//! builds standalone (see the `bsyn` feature in Cargo.toml).

/// Bitmask encoding of an IUPAC code: bit0=A, bit1=C, bit2=G, bit3=T.
/// `N` is 0b1111, an unrecognised byte is 0 and therefore matches nothing.
#[inline]
pub const fn iupac_mask(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0b0001,
        b'C' | b'c' => 0b0010,
        b'G' | b'g' => 0b0100,
        b'T' | b't' | b'U' | b'u' => 0b1000,
        b'R' | b'r' => 0b0101, // A/G
        b'Y' | b'y' => 0b1010, // C/T
        b'S' | b's' => 0b0110, // C/G
        b'W' | b'w' => 0b1001, // A/T
        b'K' | b'k' => 0b1100, // G/T
        b'M' | b'm' => 0b0011, // A/C
        b'B' | b'b' => 0b1110, // C/G/T
        b'D' | b'd' => 0b1101, // A/G/T
        b'H' | b'h' => 0b1011, // A/C/T
        b'V' | b'v' => 0b0111, // A/C/G
        b'N' | b'n' => 0b1111,
        _ => 0,
    }
}

/// True when a concrete base `base` is compatible with an IUPAC `code`.
#[inline]
pub fn iupac_matches(code: u8, base: u8) -> bool {
    let bm = iupac_mask(base);
    // A concrete base has exactly one bit set; ambiguous input never matches.
    bm.count_ones() == 1 && (iupac_mask(code) & bm) != 0
}

/// Complement of an IUPAC code, preserving degeneracy (`Y` -> `R`, `V` -> `B`, ...).
#[inline]
pub const fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' | b'U' | b'u' => b'A',
        b'R' | b'r' => b'Y',
        b'Y' | b'y' => b'R',
        b'S' | b's' => b'S',
        b'W' | b'w' => b'W',
        b'K' | b'k' => b'M',
        b'M' | b'm' => b'K',
        b'B' | b'b' => b'V',
        b'V' | b'v' => b'B',
        b'D' | b'd' => b'H',
        b'H' | b'h' => b'D',
        b'N' | b'n' => b'N',
        other => other,
    }
}

/// Reverse complement of an IUPAC string.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// A pattern is palindromic when it equals its own reverse complement; such
/// enzymes hit the same physical locus from both strands and must be deduplicated
/// during digestion (report §4.1: "回文酶的双链同位点命中已去重").
pub fn is_palindromic(pattern: &[u8]) -> bool {
    let rc = revcomp(pattern);
    rc.eq_ignore_ascii_case(pattern)
}

/// Uppercase a sequence in place, mapping every non-ACGT byte to `N`.
/// Soft-masked (lowercase) reference regions are therefore kept, not dropped —
/// masking is a repeat annotation, not a sequence-quality statement.
pub fn normalize_in_place(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        *b = match *b {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' | b'U' | b'u' => b'T',
            _ => b'N',
        };
    }
}

/// 2-bit code for a concrete base, `None` for anything ambiguous.
#[inline]
pub const fn two_bit(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Pack up to 32 bases into a `u64` (2 bits per base, 5'->3', first base in the
/// most significant used bits). Returns `None` if the tag is longer than 32 bases
/// or contains an ambiguous base.
///
/// Tags in the 16-enzyme panel are 25–33 bp; the 33 bp CspCI tag does not fit and
/// falls back to [`hash_tag`]. Both paths are collision-resistant enough for an
/// anchor table because a candidate hit is always re-verified against the stored
/// tag sequence during counting.
pub fn pack_2bit(tag: &[u8]) -> Option<u64> {
    if tag.len() > 32 {
        return None;
    }
    let mut v: u64 = 0;
    for &b in tag {
        v = (v << 2) | two_bit(b)? as u64;
    }
    Some(v)
}

/// Canonical 64-bit hash of a tag: 2-bit packing when it fits, otherwise a
/// FNV-1a fold. Length is mixed in so that tags of different enzymes with a
/// shared prefix cannot alias.
pub fn hash_tag(tag: &[u8]) -> u64 {
    let base = match pack_2bit(tag) {
        Some(v) => v,
        None => {
            let mut h: u64 = 0xcbf2_9ce4_8422_2325;
            for &b in tag {
                h ^= b as u64;
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
            h
        }
    };
    // splitmix64 finaliser: decorrelates the low bits so a plain modulo bucketing
    // (nohash-style tables) stays uniform.
    let mut z = base
        .wrapping_add((tag.len() as u64) << 58)
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Strand-canonical hash: the smaller of the forward and reverse-complement hash.
/// Used by the counter so a read matches an anchor regardless of which strand it
/// came off.
pub fn canonical_hash(tag: &[u8]) -> u64 {
    let f = hash_tag(tag);
    let r = hash_tag(&revcomp(tag));
    f.min(r)
}

/// GC fraction over concrete bases; `N` runs are excluded from the denominator.
/// Returns `None` when the window has no concrete base at all.
pub fn gc_fraction(seq: &[u8]) -> Option<f64> {
    let mut gc = 0usize;
    let mut acgt = 0usize;
    for &b in seq {
        match b {
            b'G' | b'g' | b'C' | b'c' => {
                gc += 1;
                acgt += 1;
            }
            b'A' | b'a' | b'T' | b't' | b'U' | b'u' => acgt += 1,
            _ => {}
        }
    }
    if acgt == 0 {
        None
    } else {
        Some(gc as f64 / acgt as f64)
    }
}

/// Quantise a GC fraction into the `local_gc: u8` field of an [`crate::anchor_db::Anchor`]:
/// 0..=200 encodes 0.0..=1.0 in 0.5 % steps, 255 is the "undefined" sentinel.
#[inline]
pub fn quantize_gc(frac: Option<f64>) -> u8 {
    match frac {
        Some(f) if f.is_finite() => (f.clamp(0.0, 1.0) * 200.0).round() as u8,
        _ => GC_UNDEFINED,
    }
}

/// Sentinel stored in `local_gc` when the ±flank window was all `N`.
pub const GC_UNDEFINED: u8 = 255;

/// Inverse of [`quantize_gc`].
#[inline]
pub fn dequantize_gc(q: u8) -> Option<f64> {
    if q == GC_UNDEFINED {
        None
    } else {
        Some(q as f64 / 200.0)
    }
}

/// Hamming distance, short-circuiting once `budget` is exceeded. Returns `None`
/// when the distance is greater than `budget` (or lengths differ).
///
/// This is the mismatch-budget primitive behind Syn2bANI's `tag_matcher`.
#[inline]
pub fn hamming_within(a: &[u8], b: &[u8], budget: u32) -> Option<u32> {
    if a.len() != b.len() {
        return None;
    }
    let mut d = 0u32;
    for i in 0..a.len() {
        if a[i] != b[i] {
            d += 1;
            if d > budget {
                return None;
            }
        }
    }
    Some(d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_degenerate_codes() {
        assert!(iupac_matches(b'Y', b'C'));
        assert!(iupac_matches(b'Y', b'T'));
        assert!(!iupac_matches(b'Y', b'A'));
        assert!(iupac_matches(b'V', b'G'));
        assert!(!iupac_matches(b'V', b'T'));
        assert!(iupac_matches(b'N', b'A'));
        // An ambiguous *subject* base never matches: reference `N` runs must not
        // silently create anchors.
        assert!(!iupac_matches(b'N', b'N'));
    }

    #[test]
    fn revcomp_preserves_degeneracy() {
        assert_eq!(revcomp(b"GAYNNNNNVTC"), b"GABNNNNNRTC".to_vec());
        assert_eq!(revcomp(b"ACGT"), b"ACGT".to_vec());
    }

    #[test]
    fn palindrome_detection() {
        assert!(is_palindromic(b"ACGT"));
        assert!(is_palindromic(b"GANNNNTC"));
        assert!(!is_palindromic(b"CGANNNNNNTGC"));
    }

    #[test]
    fn hash_is_strand_canonical() {
        assert_eq!(canonical_hash(b"ACGTACGT"), canonical_hash(b"ACGTACGT"));
        assert_eq!(
            canonical_hash(b"AACCGGTT"),
            canonical_hash(&revcomp(b"AACCGGTT"))
        );
        assert_ne!(canonical_hash(b"AAAAAAAA"), canonical_hash(b"AAAAAAAC"));
    }

    #[test]
    fn hash_mixes_length() {
        // "AA" and "AAAA" both pack to 0 without the length term.
        assert_ne!(hash_tag(b"AA"), hash_tag(b"AAAA"));
    }

    #[test]
    fn gc_ignores_n_runs() {
        assert_eq!(gc_fraction(b"GCGCNNNN"), Some(1.0));
        assert_eq!(gc_fraction(b"NNNN"), None);
        let q = quantize_gc(gc_fraction(b"GCAT"));
        assert_eq!(dequantize_gc(q), Some(0.5));
    }

    #[test]
    fn hamming_budget_short_circuits() {
        assert_eq!(hamming_within(b"AAAA", b"AAAA", 2), Some(0));
        assert_eq!(hamming_within(b"AAAA", b"AACA", 2), Some(1));
        assert_eq!(hamming_within(b"AAAA", b"CCCA", 2), None); // 3 > budget
        assert_eq!(hamming_within(b"AAAA", b"AAA", 2), None);
    }
}
