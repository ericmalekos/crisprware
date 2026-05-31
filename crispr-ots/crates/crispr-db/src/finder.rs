//! Site finder: walks a contig and emits PAM-flanked protospacer sites.
//!
//! Both strands are scanned. Reverse-strand sites are *re-oriented* into
//! the canonical 5'-protospacer-first form before encoding, so a single
//! `compare_mask` works for both strands at scan time. The recorded
//! `Position.offset` always points at the leftmost base on the forward
//! strand — same convention as SAM/BED.
//!
//! Windows that contain any non-ACGT byte (Ns, IUPAC ambiguity codes
//! in the genome) are skipped — FlashFry does the same via its
//! `[ACGTacgt]{20}` regex character class.

use crispr_encoding::{Enzyme, PamSide, Site};

use crate::position::{Position, Strand};

/// Plug-in receiver for sites the finder discovers. Lets the caller decide
/// where sites go (a `Vec`, a `BinTable`, a counter, an on-disk writer…).
pub trait SiteSink {
    fn accept(&mut self, site: Site, position: Position);
}

impl SiteSink for Vec<(Site, Position)> {
    fn accept(&mut self, site: Site, position: Position) {
        self.push((site, position));
    }
}

#[derive(Debug)]
pub struct SiteFinder {
    enzyme: Enzyme,
    /// Canonical PAMs in ASCII form, ready to compare byte-by-byte against
    /// a forward-strand window.
    fwd_pam_bytes: Vec<Vec<u8>>,
    /// Reverse complements of the canonical PAMs in ASCII form. Used to
    /// detect reverse-strand sites without copying the whole contig.
    rev_pam_bytes: Vec<Vec<u8>>,
}

impl SiteFinder {
    #[must_use]
    pub fn new(enzyme: Enzyme) -> Self {
        let fwd_pam_bytes: Vec<Vec<u8>> = enzyme
            .canonical_pams
            .iter()
            .map(|p| p.iter().map(|b| b.to_ascii()).collect())
            .collect();
        let rev_pam_bytes: Vec<Vec<u8>> = enzyme
            .canonical_pams
            .iter()
            .map(|p| p.iter().rev().map(|b| b.complement().to_ascii()).collect())
            .collect();
        Self {
            enzyme,
            fwd_pam_bytes,
            rev_pam_bytes,
        }
    }

    #[must_use]
    pub fn enzyme(&self) -> &Enzyme {
        &self.enzyme
    }

    /// Scan both strands of `seq`, emitting every PAM-flanked site to `sink`.
    /// `contig_id` is recorded verbatim in each emitted `Position`.
    pub fn scan(&self, seq: &[u8], contig_id: u32, sink: &mut impl SiteSink) {
        self.scan_forward(seq, contig_id, sink);
        self.scan_reverse(seq, contig_id, sink);
    }

    fn scan_forward(&self, seq: &[u8], contig_id: u32, sink: &mut impl SiteSink) {
        let pam_len = usize::from(self.enzyme.pam_len);
        let proto_len = usize::from(self.enzyme.protospacer_len);
        let total_len = pam_len + proto_len;
        if seq.len() < total_len {
            return;
        }
        let last_start = seq.len() - total_len;
        let pam_offset = match self.enzyme.pam_side {
            PamSide::ThreePrime => proto_len, // PAM after protospacer
            PamSide::FivePrime => 0,          // PAM before protospacer
        };
        for i in 0..=last_start {
            let window = &seq[i..i + total_len];
            if !is_pure_acgt(window) {
                continue;
            }
            let candidate_pam = &window[pam_offset..pam_offset + pam_len];
            if !self.fwd_pam_bytes.iter().any(|p| p == candidate_pam) {
                continue;
            }
            let site = Site::encode_ascii(window).with_count(1);
            let position = Position {
                contig_id,
                offset: u32::try_from(i).expect("offset fits in u32"),
                strand: Strand::Forward,
            };
            sink.accept(site, position);
        }
    }

    fn scan_reverse(&self, seq: &[u8], contig_id: u32, sink: &mut impl SiteSink) {
        let pam_len = usize::from(self.enzyme.pam_len);
        let proto_len = usize::from(self.enzyme.protospacer_len);
        let total_len = pam_len + proto_len;
        if seq.len() < total_len {
            return;
        }
        let last_start = seq.len() - total_len;
        // On the forward strand we look for the reverse complement of the
        // canonical PAM. For 3' PAM enzymes, rc(PAM) appears at the
        // *start* of the forward-strand window; for 5' PAM enzymes, at the
        // *end*.
        let rev_pam_offset = match self.enzyme.pam_side {
            PamSide::ThreePrime => 0,
            PamSide::FivePrime => proto_len,
        };
        for i in 0..=last_start {
            let window = &seq[i..i + total_len];
            if !is_pure_acgt(window) {
                continue;
            }
            let candidate_pam = &window[rev_pam_offset..rev_pam_offset + pam_len];
            if !self.rev_pam_bytes.iter().any(|p| p == candidate_pam) {
                continue;
            }
            // Re-orient into canonical form (PAM at the configured end, in
            // the same order as a forward-strand site would have it).
            let rc = reverse_complement(window);
            let site = Site::encode_ascii(&rc).with_count(1);
            let position = Position {
                contig_id,
                offset: u32::try_from(i).expect("offset fits in u32"),
                strand: Strand::Reverse,
            };
            sink.accept(site, position);
        }
    }
}

#[inline]
fn is_pure_acgt(window: &[u8]) -> bool {
    window
        .iter()
        .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => unreachable!(
                "reverse_complement called on non-ACGT byte (caller should pre-filter)"
            ),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crispr_encoding::Enzyme;

    fn collect_sites(finder: &SiteFinder, seq: &[u8]) -> Vec<(Site, Position)> {
        let mut sink: Vec<(Site, Position)> = Vec::new();
        finder.scan(seq, 0, &mut sink);
        sink
    }

    #[test]
    fn cas9_finds_single_forward_site() {
        let finder = SiteFinder::new(Enzyme::spcas9_ngg());
        // 20-mer protospacer ACGTACGTACGTACGTACGT, PAM AGG, total 23.
        let seq = b"ACGTACGTACGTACGTACGTAGG";
        let sites = collect_sites(&finder, seq);
        let fwd: Vec<_> = sites
            .iter()
            .filter(|(_, p)| p.strand == Strand::Forward)
            .collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].1.offset, 0);
    }

    #[test]
    fn cas9_finds_single_reverse_site() {
        let finder = SiteFinder::new(Enzyme::spcas9_ngg());
        // rev-comp of the forward case above. On forward strand we should
        // see CCT (= rc(AGG)) at the start, then 20 ACGT bases.
        let seq = b"CCTACGTACGTACGTACGTACGT";
        let sites = collect_sites(&finder, seq);
        let rev: Vec<_> = sites
            .iter()
            .filter(|(_, p)| p.strand == Strand::Reverse)
            .collect();
        assert_eq!(rev.len(), 1);
        assert_eq!(rev[0].1.offset, 0);

        // The encoded site should be the canonical forward orientation:
        // [protospacer][PAM] = ACGTACGTACGTACGTACGT + AGG.
        let decoded = rev[0].0.decode_ascii(23);
        assert_eq!(decoded, "ACGTACGTACGTACGTACGTAGG");
    }

    #[test]
    fn cas9_finds_both_strands_in_palindromic_input() {
        let finder = SiteFinder::new(Enzyme::spcas9_ngg());
        // 23-char window that happens to be a valid Cas9 site on both
        // strands: pick one with NGG at the end and CC at the start.
        //                 CC + 18 + NGG → 23 chars, but we need 20+NGG
        // Forward: ACGTACGTACGTACGTACGT AGG (site at offset 0)
        // Reverse: that whole string reverse-complemented is CCTACGTACGTACGTACGTACGT
        // — different sequence. Use a sequence that contains both patterns.
        let seq = b"ACGTACGTACGTACGTACGTAGGCCTACGTACGTACGTACGTACGT";
        let sites = collect_sites(&finder, seq);
        assert!(sites.iter().any(|(_, p)| p.strand == Strand::Forward));
        assert!(sites.iter().any(|(_, p)| p.strand == Strand::Reverse));
    }

    #[test]
    fn cas9_skips_sites_with_n_in_protospacer() {
        let finder = SiteFinder::new(Enzyme::spcas9_ngg());
        // N in the protospacer at position 10.
        let seq = b"ACGTACGTACNTACGTACGTAGG";
        let sites = collect_sites(&finder, seq);
        assert!(sites.is_empty());
    }

    #[test]
    fn cpf1_finds_forward_site() {
        let finder = SiteFinder::new(Enzyme::cpf1_tttn());
        // TTT + N + 23-mer = 27 chars.
        let seq = b"TTTACGTACGTACGTACGTACGTACGT";
        let sites = collect_sites(&finder, seq);
        let fwd: Vec<_> = sites
            .iter()
            .filter(|(_, p)| p.strand == Strand::Forward)
            .collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].1.offset, 0);
    }

    #[test]
    fn cpf1_finds_reverse_site() {
        let finder = SiteFinder::new(Enzyme::cpf1_tttn());
        // Reverse-strand Cpf1 site: forward sees rc(TTTN + 23mer)
        // = rc(23mer) + rc(TTTN). Use TTTA → TAAA at the end of the
        // forward read. 23-mer of all A's → rc = 23-mer of all T's...
        // simpler: pick a palindromic 23-mer. ACGTACGTACGTACGTACGTACG
        // → rc = CGTACGTACGTACGTACGTACGT (not palindromic). Just use a
        // hand-picked sequence and verify by decoding.
        let seq = b"ACGTACGTACGTACGTACGTACGTAAA";
        let sites = collect_sites(&finder, seq);
        let rev: Vec<_> = sites
            .iter()
            .filter(|(_, p)| p.strand == Strand::Reverse)
            .collect();
        assert_eq!(rev.len(), 1);
        // Decoded canonical form should start with TTT (the PAM) and
        // continue with the protospacer.
        let decoded = rev[0].0.decode_ascii(27);
        assert_eq!(&decoded[..3], "TTT");
    }

    #[test]
    fn short_sequence_yields_no_sites() {
        let finder = SiteFinder::new(Enzyme::spcas9_ngg());
        let sites = collect_sites(&finder, b"ACGTACGT");
        assert!(sites.is_empty());
    }

    #[test]
    fn cas9_finds_two_overlapping_forward_sites() {
        let finder = SiteFinder::new(Enzyme::spcas9_ngg());
        // Two adjacent NGGs: protospacer ends at positions 22 and 23.
        // ACGTACGTACGTACGTACGTAGG (site at 0)
        //  CGTACGTACGTACGTACGTAGGG → starting at offset 1, PAM would be at
        // pos 21-23 = "GGG", which is a valid NGG (specifically GGG). ✓
        let seq = b"ACGTACGTACGTACGTACGTAGGG";
        let sites = collect_sites(&finder, seq);
        let fwd: Vec<u32> = sites
            .iter()
            .filter(|(_, p)| p.strand == Strand::Forward)
            .map(|(_, p)| p.offset)
            .collect();
        assert_eq!(fwd, vec![0, 1]);
    }
}
