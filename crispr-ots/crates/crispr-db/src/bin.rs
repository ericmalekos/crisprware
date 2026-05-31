//! In-memory bin storage for discovered sites.
//!
//! Sites are bucketed by the first `bin_width` bases of their *protospacer*
//! (FlashFry's discovery convention, see `crispr/BinWriter.scala:59-66`).
//! For 3' PAM enzymes this is the first 7 bases of the encoded site; for
//! 5' PAM enzymes it's bases `pam_len..pam_len+7`.
//!
//! Phase 1b uses a `HashMap<u32, Vec<...>>` for clarity. Phase 1c will
//! likely swap in a dense `Vec` indexed by bin id — the keyspace is small
//! (`4^7 = 16_384` bins for a 7-mer prefix), and string-keyed hashing
//! repeated per bin access was a perf hazard we explicitly called out
//! in `NOTES.md`.

use std::collections::HashMap;
use std::io;
use std::path::Path;

use crispr_encoding::{Enzyme, PamSide, Site};
use serde::{Deserialize, Serialize};

use crate::finder::SiteSink;
use crate::mmap_db::write_binary_db;
use crate::packed::PackedEntry;
use crate::position::Position;

const DEFAULT_BIN_WIDTH: u8 = 7;

#[derive(Debug, Serialize, Deserialize)]
pub struct BinTable {
    enzyme: Enzyme,
    bin_width: u8,
    /// Contig names registered during build, indexed by the `contig_id`
    /// stored in each [`Position`]. Used by output layers to recover the
    /// FASTA-header name for a given site.
    contig_names: Vec<String>,
    /// Per-bin lists of packed `(site, position)` records. The
    /// `PackedEntry` representation is shared with the v2 on-disk format
    /// so the in-memory and mmap'd paths can use identical Scanner code.
    bins: HashMap<u32, Vec<PackedEntry>>,
}

impl BinTable {
    /// Build a bin table for `enzyme` with the default 7-bp bin width.
    /// Takes the `Enzyme` by value — callers that share the enzyme with a
    /// `SiteFinder` should `.clone()` it.
    #[must_use]
    pub fn for_enzyme(enzyme: Enzyme) -> Self {
        Self::with_width(enzyme, DEFAULT_BIN_WIDTH)
    }

    /// Build a bin table for `enzyme` with a custom bin width.
    ///
    /// # Panics
    /// Panics if `bin_width` is 0, exceeds the enzyme's protospacer length,
    /// or exceeds 15 (the u32 key limit, beyond which `4^bin_width` no
    /// longer fits in a `u32`).
    #[must_use]
    pub fn with_width(enzyme: Enzyme, bin_width: u8) -> Self {
        assert!(
            bin_width <= enzyme.protospacer_len,
            "bin width {bin_width} exceeds protospacer length {}",
            enzyme.protospacer_len
        );
        assert!(bin_width > 0, "bin width must be at least 1");
        assert!(
            bin_width <= 15,
            "bin width {bin_width} exceeds the 15-base u32 key limit"
        );
        Self {
            enzyme,
            bin_width,
            contig_names: Vec::new(),
            bins: HashMap::new(),
        }
    }

    /// Reference to the enzyme this table was built for. Scanners pull
    /// `compare_mask`, `protospacer_len`, etc. from here.
    #[must_use]
    pub fn enzyme(&self) -> &Enzyme {
        &self.enzyme
    }

    /// Bin prefix width in bases.
    #[must_use]
    pub fn bin_width(&self) -> u8 {
        self.bin_width
    }

    /// Register a contig name and return the `u32` id to pass to
    /// [`crate::SiteFinder::scan`]. Ids are dense (0, 1, 2, …) in the
    /// order contigs are registered.
    ///
    /// # Panics
    /// Panics if the contig count exceeds `u32::MAX`.
    pub fn register_contig(&mut self, name: impl Into<String>) -> u32 {
        let id = u32::try_from(self.contig_names.len()).expect("contig count fits in u32");
        self.contig_names.push(name.into());
        id
    }

    /// Look up a contig name by id, returning `None` if the id was never
    /// registered (e.g. tests that scan with an ad-hoc id of 0).
    #[must_use]
    pub fn contig_name(&self, id: u32) -> Option<&str> {
        self.contig_names.get(id as usize).map(String::as_str)
    }

    /// Number of registered contigs.
    #[must_use]
    pub fn contig_count(&self) -> usize {
        self.contig_names.len()
    }

    /// Serialize the table to `path` in the crispr-ots v2 database format
    /// (header + per-bin offset table + packed entries; mmap-friendly).
    ///
    /// The reader for this format is [`crate::MmapDb`]. There is no
    /// in-memory `load_from_path` — re-reading a v2 file always goes
    /// through `MmapDb::open`, which mmaps the file and returns
    /// zero-copy views into the page cache rather than allocating a
    /// fresh `BinTable`. Tests that need an in-memory `BinTable` should
    /// build it via the `crispr-cli::build_table_in_memory` helper.
    ///
    /// # Errors
    /// Surfaces filesystem and serialization errors as `io::Error`.
    pub fn save_to_path(&self, path: &Path) -> io::Result<()> {
        write_binary_db(self, path)
    }

    /// Total number of sites currently stored.
    #[must_use]
    pub fn total(&self) -> usize {
        self.bins.values().map(Vec::len).sum()
    }

    /// Number of bin keys that have at least one site.
    #[must_use]
    pub fn occupied_bins(&self) -> usize {
        self.bins.len()
    }

    /// Iterate `(bin_key, packed_entries)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = (u32, &[PackedEntry])> + '_ {
        self.bins.iter().map(|(&k, v)| (k, v.as_slice()))
    }

    /// Inspect a single bin.
    #[must_use]
    pub fn get(&self, key: u32) -> Option<&[PackedEntry]> {
        self.bins.get(&key).map(Vec::as_slice)
    }

    /// Compute the bin key for `site` under this table's enzyme. Public so
    /// callers can ask "which bin would this site go into?" without
    /// actually inserting it (e.g. for inverse lookups during query).
    ///
    /// # Panics
    /// Panics if the computed key doesn't fit in a `u32` — this is
    /// statically impossible given the `bin_width <= 15` check enforced at
    /// `BinTable` construction, but the `try_from` keeps the proof local.
    #[must_use]
    pub fn bin_key(&self, site: Site) -> u32 {
        // The first protospacer base sits at the highest used bit positions:
        //   - 3' PAM (Cas9): encoded order is [protospacer][PAM]; the
        //     protospacer's first base is at the top of the encoding.
        //   - 5' PAM (Cpf1): encoded order is [PAM][protospacer]; the
        //     protospacer's first base sits *below* the PAM at the top.
        //
        // In both cases the bin (the first `bin_width` protospacer bases)
        // ends at a bit position determined by `protospacer_len` for 5'
        // PAM and by `protospacer_len + pam_len` for 3' PAM, and starts
        // `bin_width * 2` bits below that.
        let bin_bits = u32::from(self.bin_width) * 2;
        let bin_top_from_lsb = match self.enzyme.pam_side {
            PamSide::ThreePrime => {
                (u32::from(self.enzyme.protospacer_len) + u32::from(self.enzyme.pam_len)) * 2
            }
            PamSide::FivePrime => u32::from(self.enzyme.protospacer_len) * 2,
        };
        let start_bit = bin_top_from_lsb - bin_bits;
        let extracted = extract_bits(site, start_bit, bin_bits);
        u32::try_from(extracted).expect("bin key fits in u32 (bin_width <= 15)")
    }
}

impl SiteSink for BinTable {
    fn accept(&mut self, site: Site, position: Position) {
        let key = self.bin_key(site);
        self.bins
            .entry(key)
            .or_default()
            .push(PackedEntry::from_site_position(site, position));
    }
}

/// Extract `num_bits` bits starting at `start_bit` from the composite
/// (low + high) 128-bit `Site` encoding. Handles the case where the
/// requested range straddles the 48-bit low/high boundary.
fn extract_bits(site: Site, start_bit: u32, num_bits: u32) -> u64 {
    const LOW_WIDTH: u32 = 48;
    let end_bit = start_bit + num_bits;
    if end_bit <= LOW_WIDTH {
        (site.low >> start_bit) & ((1u64 << num_bits) - 1)
    } else if start_bit >= LOW_WIDTH {
        // Entirely inside `high`. `high` bits 0-7 hold the sequence; we
        // don't shift any further to avoid masking the count field
        // (callers should only request sequence bits).
        let local_shift = start_bit - LOW_WIDTH;
        (site.high >> local_shift) & ((1u64 << num_bits) - 1)
    } else {
        // Straddle: lower portion from `low`, upper portion from `high`.
        let low_count = LOW_WIDTH - start_bit;
        let high_count = num_bits - low_count;
        let low_part = site.low >> start_bit;
        let high_part = (site.high & ((1u64 << high_count) - 1)) << low_count;
        low_part | high_part
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::position::Strand;
    use crispr_encoding::{Base, Enzyme};

    /// Build a u32 bin key by packing `bases` in the same FlashFry "first
    /// base at the highest used bits" convention used by `Site` encoding.
    fn pack_bin_key(bases: &[Base]) -> u32 {
        let mut k: u32 = 0;
        for b in bases {
            k = (k << 2) | u32::from(b.bits());
        }
        k
    }

    #[test]
    fn cas9_bin_key_is_first_7_protospacer_bases() {
        let enz = Enzyme::spcas9_ngg();
        let table = BinTable::for_enzyme(enz);
        // 20-mer protospacer = ACGTACGAAAAAAAAAAAAAA, PAM AGG.
        let site = Site::encode_ascii(b"ACGTACGAAAAAAAAAAAAAAGG");
        let key = table.bin_key(site);
        let expected = pack_bin_key(&[
            Base::A,
            Base::C,
            Base::G,
            Base::T,
            Base::A,
            Base::C,
            Base::G,
        ]);
        assert_eq!(key, expected);
    }

    #[test]
    fn cpf1_bin_key_skips_pam() {
        let enz = Enzyme::cpf1_tttn();
        let table = BinTable::for_enzyme(enz);
        // TTT + N=A + 20-mer protospacer ACGTACGAAAAAAAAAAAAA.
        let site = Site::encode_ascii(b"TTTAACGTACGAAAAAAAAAAAAA");
        assert_eq!(site.decode_ascii(24).len(), 24);
        let key = table.bin_key(site);
        // Bin key should be the first 7 protospacer bases ACGTACG,
        // *not* including any of TTT[N].
        let expected = pack_bin_key(&[
            Base::A,
            Base::C,
            Base::G,
            Base::T,
            Base::A,
            Base::C,
            Base::G,
        ]);
        assert_eq!(key, expected);
    }

    #[test]
    fn bin_key_straddles_low_high_for_long_3prime_pam_enzyme() {
        // Synthetic enzyme: 3' PAM, 23-nt protospacer, 3-nt PAM, total 26.
        // The bin (first 7 protospacer bases) sits at the top of the
        // encoding: starts at bit (26 - 7) * 2 = 38, ends at bit 52.
        // `low` holds bits 0-47 — so bits 38-47 (5 bases) come from low,
        // bits 48-51 (2 bases) come from high. Exercises the straddle path.
        let enz = Enzyme::from_iupac(
            "TestLongCas9",
            &[
                crispr_encoding::IupacCode::N,
                crispr_encoding::IupacCode::G,
                crispr_encoding::IupacCode::G,
            ],
            PamSide::ThreePrime,
            23,
        );
        let table = BinTable::for_enzyme(enz);
        // 23-base protospacer (ACGTACG + 16 A's) + AGG PAM.
        // 7 + 16 + 3 = 26 chars. The PAM's leading A happens to look like
        // a 17th protospacer A when read inline, so eyeball-counting the A's
        // is misleading — trust the constants.
        let seq = b"ACGTACGAAAAAAAAAAAAAAAAAGG"; // 7 + 17 + 2 = 26
        assert_eq!(seq.len(), 26);
        let site = Site::encode_ascii(seq);
        let key = table.bin_key(site);
        let expected = pack_bin_key(&[
            Base::A,
            Base::C,
            Base::G,
            Base::T,
            Base::A,
            Base::C,
            Base::G,
        ]);
        assert_eq!(key, expected);
    }

    #[test]
    fn bin_table_groups_sites_by_key() {
        let enz = Enzyme::spcas9_ngg();
        let mut table = BinTable::for_enzyme(enz);
        let pos_a = Position {
            contig_id: 0,
            offset: 0,
            strand: Strand::Forward,
        };
        let pos_b = Position {
            contig_id: 0,
            offset: 100,
            strand: Strand::Forward,
        };
        let pos_c = Position {
            contig_id: 0,
            offset: 200,
            strand: Strand::Forward,
        };
        // Two sites that share the same 7-mer prefix ACGTACG, one with a
        // different prefix.
        let site_a = Site::encode_ascii(b"ACGTACGAAAAAAAAAAAAAAGG");
        let site_b = Site::encode_ascii(b"ACGTACGCCCCCCCCCCCCCAGG");
        let site_c = Site::encode_ascii(b"TTTTTTTAAAAAAAAAAAAAAGG");
        table.accept(site_a, pos_a);
        table.accept(site_b, pos_b);
        table.accept(site_c, pos_c);
        assert_eq!(table.total(), 3);
        assert_eq!(table.occupied_bins(), 2);

        let shared_key = table.bin_key(site_a);
        assert_eq!(table.bin_key(site_b), shared_key);
        assert_eq!(table.get(shared_key).expect("bin exists").len(), 2);
    }

    #[test]
    fn cpf1_bin_key_independent_of_pam_n() {
        // Two sites with the same protospacer but different PAM N — the
        // bin key (first 7 protospacer bases) should be identical because
        // the PAM lives outside the bin window.
        let enz = Enzyme::cpf1_tttn();
        let table = BinTable::for_enzyme(enz);
        let site_x = Site::encode_ascii(b"TTTAACGTACGAAAAAAAAAAAAA");
        let site_y = Site::encode_ascii(b"TTTCACGTACGAAAAAAAAAAAAA");
        assert_eq!(table.bin_key(site_x), table.bin_key(site_y));
    }

    #[test]
    #[should_panic(expected = "bin width 30 exceeds protospacer length 20")]
    fn bin_width_overflow_panics() {
        let _ = BinTable::with_width(Enzyme::spcas9_ngg(), 30);
    }

    #[test]
    fn persistence_roundtrip() {
        // BinTable -> v2 on-disk -> MmapDb is exercised end-to-end in
        // `crate::mmap_db::tests::roundtrip_basic`. Keep this stub as a
        // smoke check that the BinTable side of save_to_path doesn't
        // panic on the typical input shape.
        use crate::MmapDb;
        let dir = std::env::temp_dir().join(format!(
            "crispr-db-persist-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_nanos()),
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.crot2");

        let mut table = BinTable::for_enzyme(Enzyme::spcas9_ngg());
        table.register_contig("chrTest".to_string());
        for (seq, off, st) in [
            (b"ACGTACGAAAAAAAAAAAAAAGG", 10_u32, Strand::Forward),
            (b"ACGTACGCCCCCCCCCCCCCAGG", 20, Strand::Forward),
            (b"TTTTTTTTAAAAAAAAAAAAAGG", 30, Strand::Reverse),
        ] {
            table.accept(
                Site::encode_ascii(seq),
                Position {
                    contig_id: 0,
                    offset: off,
                    strand: st,
                },
            );
        }

        table.save_to_path(&path).expect("save");
        let mm = MmapDb::open(&path).expect("open");
        assert_eq!(mm.contig_name(0), Some("chrTest"));
        let total: usize = mm.iter().map(|(_, e)| e.len()).sum();
        assert_eq!(total, table.total());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
