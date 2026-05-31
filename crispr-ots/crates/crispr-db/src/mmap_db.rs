// The reader/writer freely converts between u64 (on-disk widths) and
// usize (in-memory indexing). On any 64-bit target — which is all we
// support — these are width-equivalent; on hypothetical 32-bit targets
// the cast would truncate, but the runtime bounds checks in `open()`
// catch any out-of-range value before it's used.
#![allow(clippy::cast_possible_truncation)]

//! Zero-copy reader (and writer) for the v2 on-disk database format.
//!
//! # File layout
//!
//! ```text
//! +--------------------------------+ offset 0
//! | DbHeader (128 bytes, repr(C))  |
//! +--------------------------------+ enzyme_name_off
//! | enzyme name (UTF-8, raw)       |
//! +--------------------------------+ pam_data_off
//! | canonical PAMs (pam_count ×    |
//! |   pam_len bytes, ASCII)        |
//! +--------------------------------+ contig_table_off
//! | contig table:                  |
//! |   for each contig:             |
//! |     u32 name_len               |
//! |     name bytes (UTF-8)         |
//! +--------------------------------+ bin_offsets_off (8-byte aligned)
//! | BinOffset[bin_count]           |
//! |   (start_index: u32,           |
//! |    count: u32)                 |
//! +--------------------------------+ entries_off (8-byte aligned)
//! | PackedEntry[total_entries]     |
//! +--------------------------------+ EOF
//! ```
//!
//! The bin offsets table is indexed directly by `bin_key`. Empty bins
//! have `count == 0`. `start_index` is a *count* into the entries array,
//! not a byte offset, so a single bin's payload is
//! `entries[off.start_index..off.start_index + off.count]`.
//!
//! `bytemuck::Pod` lets us treat the mmap'd byte region as a typed slice
//! of `BinOffset` / `PackedEntry` with zero deserialization cost.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::mem::size_of;
use std::path::Path;

use bytemuck::{Pod, Zeroable};
use crispr_encoding::{Base, Enzyme, PamSide, SiteMask};
use memmap2::Mmap;

use crate::bin::BinTable;
use crate::packed::PackedEntry;

const MAGIC: [u8; 8] = *b"CRSPROT2";
const VERSION: u32 = 2;

/// On-disk database header. Always exactly 128 bytes, repr(C), little-
/// endian fields. The trailing `_reserved` padding keeps us at a clean
/// 128 so future fields can be added without bumping the format version
/// (callers that don't know about them just ignore the bytes).
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
#[repr(C)]
pub struct DbHeader {
    magic: [u8; 8],
    version: u32,
    flags: u32,

    bin_width: u32,
    pam_side: u32,
    pam_len: u32,
    protospacer_len: u32,

    pam_count: u32,
    contig_count: u32,
    bin_count: u32,
    _pad0: u32,

    total_entries: u64,
    compare_mask_low: u64,
    compare_mask_high: u64,

    enzyme_name_off: u64,
    pam_data_off: u64,
    contig_table_off: u64,
    bin_offsets_off: u64,
    entries_off: u64,

    // Pad the header out to exactly 128 bytes so future fields can be
    // appended without bumping the format version.
    _reserved: [u64; 2],
}

const _: () = {
    assert!(size_of::<DbHeader>() == 128);
};

/// One entry of the bin-offsets table. Both fields are counts within the
/// `entries` section, not byte offsets — the entries are all 24 bytes,
/// so `byte_off = start_index * 24`.
#[derive(Clone, Copy, Debug, Default, Pod, Zeroable)]
#[repr(C)]
pub struct BinOffset {
    pub start_index: u32,
    pub count: u32,
}

const _: () = {
    assert!(size_of::<BinOffset>() == 8);
};

/// Write a `BinTable` to disk in the v2 format.
///
/// # Errors
/// Surfaces I/O errors from the underlying file.
///
/// # Panics
/// Panics if any quantity overflows the on-disk width (e.g. more than
/// `u32::MAX` total entries, more than 1 G bins). All current enzyme/
/// genome combinations stay well within these limits.
#[allow(clippy::too_many_lines)] // single-shot layout + write reads top-to-bottom
pub fn write_binary_db(table: &BinTable, path: &Path) -> io::Result<()> {
    let enzyme = table.enzyme();
    let bin_width = table.bin_width();
    let bin_count = 4u64
        .checked_pow(u32::from(bin_width))
        .expect("bin count fits u64");

    // Collect bins in canonical (bin_key) order.
    let mut bins_in_order: Vec<(u32, &[PackedEntry])> = table.iter().collect();
    bins_in_order.sort_unstable_by_key(|&(k, _)| k);

    let total_entries: u64 = bins_in_order.iter().map(|(_, e)| e.len() as u64).sum();
    let total_entries_u32 =
        u32::try_from(total_entries).expect("total entries fits in u32 (< 4 G sites)");

    // ---- Build the bin offsets table. ----
    // Dense Vec indexed by bin_key; bins not present in the source have
    // count == 0 (and arbitrary start_index — we set it to a sane value
    // for forward iteration in case a reader uses it).
    let bin_count_us = usize::try_from(bin_count).expect("bin count fits usize");
    let mut bin_offsets: Vec<BinOffset> = vec![BinOffset::default(); bin_count_us];
    let mut cursor: u32 = 0;
    // First populate occupied bins.
    for (key, entries) in &bins_in_order {
        let len = u32::try_from(entries.len()).expect("bin entries fit u32");
        bin_offsets[*key as usize] = BinOffset {
            start_index: cursor,
            count: len,
        };
        cursor = cursor.checked_add(len).expect("cursor fits u32");
    }
    // Fill unoccupied bins with start_index = end-of-section so that a
    // forward iterator that doesn't check `count` doesn't read garbage.
    for bo in &mut bin_offsets {
        if bo.count == 0 {
            bo.start_index = total_entries_u32;
        }
    }

    // ---- Lay out the byte sections. ----
    let enzyme_name_bytes = enzyme.name.as_bytes();
    let pam_count_u32 =
        u32::try_from(enzyme.canonical_pams.len()).expect("PAM count fits u32 (typical: 1–16)");
    let pam_len_u8 = enzyme.pam_len;
    let pam_data: Vec<u8> = enzyme
        .canonical_pams
        .iter()
        .flat_map(|pam| pam.iter().map(|b| b.to_ascii()))
        .collect();
    let mut contig_table: Vec<u8> = Vec::new();
    for i in 0..table.contig_count() {
        let name = table.contig_name(u32::try_from(i).unwrap()).unwrap_or("");
        let len_bytes = u32::try_from(name.len())
            .expect("contig name length fits u32")
            .to_le_bytes();
        contig_table.extend_from_slice(&len_bytes);
        contig_table.extend_from_slice(name.as_bytes());
    }

    let mut off: u64 = size_of::<DbHeader>() as u64;
    let enzyme_name_off = off;
    off += enzyme_name_bytes.len() as u64;
    let pam_data_off = off;
    off += pam_data.len() as u64;
    let contig_table_off = off;
    off += contig_table.len() as u64;
    // 8-byte align before the BinOffset table.
    let bin_offsets_off = (off + 7) & !7;
    let pad_before_offsets = bin_offsets_off - off;
    off = bin_offsets_off + bin_count * (size_of::<BinOffset>() as u64);
    // 8-byte align before the entries array (already 8-aligned since
    // BinOffset is 8 bytes wide, but keep the round for safety).
    let entries_off = (off + 7) & !7;
    let pad_before_entries = entries_off - off;

    // ---- Build header. ----
    let pam_side_u32 = match enzyme.pam_side {
        PamSide::ThreePrime => 0_u32,
        PamSide::FivePrime => 1,
    };
    let header = DbHeader {
        magic: MAGIC,
        version: VERSION,
        flags: 0,
        bin_width: u32::from(bin_width),
        pam_side: pam_side_u32,
        pam_len: u32::from(pam_len_u8),
        protospacer_len: u32::from(enzyme.protospacer_len),
        pam_count: pam_count_u32,
        contig_count: u32::try_from(table.contig_count()).expect("contig count fits u32"),
        bin_count: u32::try_from(bin_count).expect("bin count fits u32"),
        _pad0: 0,
        total_entries,
        compare_mask_low: enzyme.compare_mask.low,
        compare_mask_high: enzyme.compare_mask.high,
        enzyme_name_off,
        pam_data_off,
        contig_table_off,
        bin_offsets_off,
        entries_off,
        _reserved: [0; 2],
    };

    // ---- Write. ----
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    w.write_all(bytemuck::bytes_of(&header))?;
    w.write_all(enzyme_name_bytes)?;
    w.write_all(&pam_data)?;
    w.write_all(&contig_table)?;
    if pad_before_offsets > 0 {
        w.write_all(&vec![0_u8; pad_before_offsets as usize])?;
    }
    w.write_all(bytemuck::cast_slice(&bin_offsets))?;
    if pad_before_entries > 0 {
        w.write_all(&vec![0_u8; pad_before_entries as usize])?;
    }
    // Write entries in canonical bin order — same ordering used to
    // populate `bin_offsets.start_index` above.
    for (_, entries) in &bins_in_order {
        w.write_all(bytemuck::cast_slice(entries))?;
    }
    w.flush()?;
    Ok(())
}

/// Memory-mapped view of a v2 database. Constructed once at discover
/// startup and queried zero-copy.
#[derive(Debug)]
pub struct MmapDb {
    // The mmap must outlive every slice we hand out; keep it as the first
    // field so it's dropped last.
    mmap: Mmap,
    enzyme: Enzyme,
    bin_width: u8,
    bin_count: u32,
    contig_names: Vec<String>,
    bin_offsets_ptr_off: usize,
    entries_ptr_off: usize,
}

/// Error variants from `MmapDb::open`.
#[derive(Debug)]
pub enum MmapDbError {
    Io(io::Error),
    BadMagic,
    UnsupportedVersion(u32),
    Truncated(&'static str),
}

impl std::fmt::Display for MmapDbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::BadMagic => write!(f, "not a crispr-ots v2 database (magic mismatch)"),
            Self::UnsupportedVersion(v) => {
                write!(f, "unsupported database version v{v} (expected v{VERSION})")
            }
            Self::Truncated(section) => write!(f, "database truncated in {section} section"),
        }
    }
}

impl std::error::Error for MmapDbError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for MmapDbError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

impl MmapDb {
    /// Open a v2 database for read-only zero-copy access.
    ///
    /// # Errors
    /// Returns [`MmapDbError`] for I/O, bad magic, version mismatch, or
    /// any structural inconsistency caught at open time.
    pub fn open(path: &Path) -> Result<Self, MmapDbError> {
        let file = File::open(path)?;
        // SAFETY: the file is opened read-only; we never write into the
        // mmap; concurrent modification of the file by another process
        // while we read is undefined per memmap2's docs but accepted as
        // the standard trade-off for zero-copy database files.
        let mmap = unsafe { Mmap::map(&file)? };
        Self::from_mmap(mmap)
    }

    fn from_mmap(mmap: Mmap) -> Result<Self, MmapDbError> {
        if mmap.len() < size_of::<DbHeader>() {
            return Err(MmapDbError::Truncated("header"));
        }
        let header: DbHeader = *bytemuck::from_bytes(&mmap[..size_of::<DbHeader>()]);
        if header.magic != MAGIC {
            return Err(MmapDbError::BadMagic);
        }
        if header.version != VERSION {
            return Err(MmapDbError::UnsupportedVersion(header.version));
        }
        if header.bin_width > 15 {
            return Err(MmapDbError::Truncated("header: bin_width out of range"));
        }

        // Bounds-check every section against mmap length up front.
        let enzyme_name_end = (header.enzyme_name_off as usize) + enzyme_name_len(&header);
        if mmap.len() < enzyme_name_end {
            return Err(MmapDbError::Truncated("enzyme name"));
        }
        let pam_data_len = (header.pam_count as usize) * (header.pam_len as usize);
        let pam_data_end = (header.pam_data_off as usize) + pam_data_len;
        if mmap.len() < pam_data_end {
            return Err(MmapDbError::Truncated("PAM data"));
        }
        let bin_offsets_bytes = (header.bin_count as usize) * size_of::<BinOffset>();
        let bin_offsets_end = (header.bin_offsets_off as usize) + bin_offsets_bytes;
        if mmap.len() < bin_offsets_end {
            return Err(MmapDbError::Truncated("bin offsets"));
        }
        let entries_bytes = (header.total_entries as usize) * size_of::<PackedEntry>();
        let entries_end = (header.entries_off as usize) + entries_bytes;
        if mmap.len() < entries_end {
            return Err(MmapDbError::Truncated("entries"));
        }

        // Parse the enzyme metadata into an owned `Enzyme`. This data
        // lives in the mmap region, but `Enzyme` predates the v2 format
        // and owns its strings/`Vec`s — copying it once at open time is
        // cheap (typically a few hundred bytes).
        let enzyme_name =
            std::str::from_utf8(&mmap[(header.enzyme_name_off as usize)..enzyme_name_end])
                .map_err(|_| MmapDbError::Truncated("enzyme name UTF-8"))?
                .to_string();

        let pam_section = &mmap[(header.pam_data_off as usize)..pam_data_end];
        let pam_chunk_size = header.pam_len as usize;
        let canonical_pams: Vec<Vec<Base>> = pam_section
            .chunks_exact(pam_chunk_size)
            .map(|chunk| {
                chunk
                    .iter()
                    .map(|&c| Base::from_ascii(c).ok_or(MmapDbError::Truncated("PAM byte")))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let pam_side = match header.pam_side {
            0 => PamSide::ThreePrime,
            1 => PamSide::FivePrime,
            _ => return Err(MmapDbError::Truncated("PAM side")),
        };
        // Recompute compare_mask from the stored fields, ignoring what's
        // in the header (it's there for sanity-check / future use).
        let pam_len_u8 =
            u8::try_from(header.pam_len).map_err(|_| MmapDbError::Truncated("pam_len"))?;
        let protospacer_len_u8 = u8::try_from(header.protospacer_len)
            .map_err(|_| MmapDbError::Truncated("protospacer_len"))?;
        let enzyme = Enzyme::from_canonical(
            enzyme_name,
            canonical_pams,
            pam_side,
            pam_len_u8,
            protospacer_len_u8,
        );
        // Cross-check the stored compare_mask matches what `from_canonical`
        // derived. Drift here would indicate a corrupted file.
        let expected_mask = SiteMask {
            low: header.compare_mask_low,
            high: header.compare_mask_high,
        };
        if enzyme.compare_mask != expected_mask {
            return Err(MmapDbError::Truncated("compare_mask drift"));
        }

        let contig_names = parse_contig_table(
            &mmap[(header.contig_table_off as usize)..(header.bin_offsets_off as usize)],
            header.contig_count as usize,
        )?;

        Ok(Self {
            mmap,
            enzyme,
            bin_width: u8::try_from(header.bin_width).expect("bin_width fits u8"),
            bin_count: header.bin_count,
            contig_names,
            bin_offsets_ptr_off: header.bin_offsets_off as usize,
            entries_ptr_off: header.entries_off as usize,
        })
    }

    /// Total number of bins (including empties).
    #[must_use]
    pub fn bin_count(&self) -> u32 {
        self.bin_count
    }

    #[must_use]
    pub fn enzyme(&self) -> &Enzyme {
        &self.enzyme
    }

    #[must_use]
    pub fn bin_width(&self) -> u8 {
        self.bin_width
    }

    #[must_use]
    pub fn contig_name(&self, id: u32) -> Option<&str> {
        self.contig_names.get(id as usize).map(String::as_str)
    }

    #[must_use]
    pub fn contig_count(&self) -> usize {
        self.contig_names.len()
    }

    /// Full bin-offsets table, indexed by `bin_key`.
    #[must_use]
    pub fn bin_offsets(&self) -> &[BinOffset] {
        let len = self.bin_count as usize;
        let bytes = &self.mmap
            [self.bin_offsets_ptr_off..self.bin_offsets_ptr_off + len * size_of::<BinOffset>()];
        bytemuck::cast_slice(bytes)
    }

    /// All entries concatenated in canonical bin order.
    #[must_use]
    pub fn entries(&self) -> &[PackedEntry] {
        bytemuck::cast_slice(&self.mmap[self.entries_ptr_off..])
    }

    /// Slice of entries for a specific bin. `None` if the bin is empty.
    #[must_use]
    pub fn bin(&self, key: u32) -> Option<&[PackedEntry]> {
        let off = *self.bin_offsets().get(key as usize)?;
        if off.count == 0 {
            return None;
        }
        let start = off.start_index as usize;
        let end = start + off.count as usize;
        Some(&self.entries()[start..end])
    }

    /// Iterate `(bin_key, entries)` pairs over non-empty bins — same
    /// shape as `BinTable::iter`.
    ///
    /// # Panics
    /// Panics if a bin's `start_index + count` overflows `usize`, which
    /// cannot happen on any 64-bit target given the `total_entries`
    /// bounds check in `open()`.
    pub fn iter(&self) -> impl Iterator<Item = (u32, &[PackedEntry])> + '_ {
        self.bin_offsets()
            .iter()
            .enumerate()
            .filter(|(_, off)| off.count > 0)
            .map(|(k, off)| {
                let start = off.start_index as usize;
                let end = start + off.count as usize;
                let entries: &[PackedEntry] = &self.entries()[start..end];
                (u32::try_from(k).expect("bin index fits u32"), entries)
            })
    }

    /// Mirror of `BinTable::bin_key`: compute which bin a guide would land
    /// in, given the database's enzyme + bin width. Used by the scanner's
    /// guide-prefix prefilter.
    #[must_use]
    pub fn bin_key(&self, site: crispr_encoding::Site) -> u32 {
        let bin_bits = u32::from(self.bin_width) * 2;
        let bin_top_from_lsb = match self.enzyme.pam_side {
            PamSide::ThreePrime => {
                (u32::from(self.enzyme.protospacer_len) + u32::from(self.enzyme.pam_len)) * 2
            }
            PamSide::FivePrime => u32::from(self.enzyme.protospacer_len) * 2,
        };
        let start_bit = bin_top_from_lsb - bin_bits;
        extract_bits(site, start_bit, bin_bits)
    }
}

fn enzyme_name_len(header: &DbHeader) -> usize {
    (header.pam_data_off - header.enzyme_name_off) as usize
}

fn parse_contig_table(bytes: &[u8], expected_count: usize) -> Result<Vec<String>, MmapDbError> {
    let mut out = Vec::with_capacity(expected_count);
    let mut cur = 0;
    while out.len() < expected_count {
        if bytes.len() < cur + 4 {
            return Err(MmapDbError::Truncated("contig table length prefix"));
        }
        let len_bytes: [u8; 4] = bytes[cur..cur + 4].try_into().unwrap();
        let len = u32::from_le_bytes(len_bytes) as usize;
        cur += 4;
        if bytes.len() < cur + len {
            return Err(MmapDbError::Truncated("contig name body"));
        }
        let name = std::str::from_utf8(&bytes[cur..cur + len])
            .map_err(|_| MmapDbError::Truncated("contig name UTF-8"))?
            .to_string();
        out.push(name);
        cur += len;
    }
    Ok(out)
}

/// Bit-extract helper, copy of the one in `bin.rs`. Inlined here so the
/// mmap reader doesn't have to depend on internal `bin` module symbols.
fn extract_bits(site: crispr_encoding::Site, start_bit: u32, num_bits: u32) -> u32 {
    const LOW_WIDTH: u32 = 48;
    let end_bit = start_bit + num_bits;
    let raw = if end_bit <= LOW_WIDTH {
        (site.low >> start_bit) & ((1u64 << num_bits) - 1)
    } else if start_bit >= LOW_WIDTH {
        let local_shift = start_bit - LOW_WIDTH;
        (site.high >> local_shift) & ((1u64 << num_bits) - 1)
    } else {
        let low_count = LOW_WIDTH - start_bit;
        let high_count = num_bits - low_count;
        let low_part = site.low >> start_bit;
        let high_part = (site.high & ((1u64 << high_count) - 1)) << low_count;
        low_part | high_part
    };
    u32::try_from(raw).expect("bin key fits in u32")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::position::{Position, Strand};
    use crispr_encoding::{Enzyme, Site};

    fn tmp_dir() -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!(
            "crispr-ots-mmap-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_nanos())
        ));
        std::fs::create_dir_all(&p).unwrap();
        p
    }

    #[test]
    fn roundtrip_basic() {
        let dir = tmp_dir();
        let path = dir.join("a.crot2");
        let mut table = BinTable::for_enzyme(Enzyme::spcas9_ngg());
        table.register_contig("chrA".to_string());
        for (seq, off, st) in [
            (b"ACGTACGAAAAAAAAAAAAAAGG", 0_u32, Strand::Forward),
            (b"ACGTACGCCCCCCCCCCCCCAGG", 50, Strand::Forward),
            (b"TTTTTTTTAAAAAAAAAAAAAGG", 100, Strand::Reverse),
            (b"ACGTACGAAAAAAAAAAAAAAGG", 200, Strand::Forward),
        ] {
            use crate::SiteSink;
            table.accept(
                Site::encode_ascii(seq),
                Position {
                    contig_id: 0,
                    offset: off,
                    strand: st,
                },
            );
        }
        write_binary_db(&table, &path).unwrap();
        let db = MmapDb::open(&path).unwrap();
        assert_eq!(db.bin_width(), table.bin_width());
        assert_eq!(db.contig_name(0), Some("chrA"));
        let total: usize = db.iter().map(|(_, e)| e.len()).sum();
        assert_eq!(total, table.total());
        // Every bin we wrote is recoverable.
        for (key, expected) in table.iter() {
            let got = db.bin(key).expect("bin present");
            let mut a = expected.to_vec();
            let mut b = got.to_vec();
            a.sort_by_key(|e| (e.position().contig_id, e.position().offset));
            b.sort_by_key(|e| (e.position().contig_id, e.position().offset));
            assert_eq!(a, b);
        }
        // bin_key on a site that matches `expected` returns the same key.
        for (key, entries) in db.iter() {
            for entry in entries {
                assert_eq!(db.bin_key(entry.site()), key);
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn open_rejects_bad_magic() {
        let dir = tmp_dir();
        let path = dir.join("bogus.crot2");
        // Need at least one header's worth of bytes or we trip the
        // truncated-header check first; pad with zeros.
        let mut bytes = vec![0_u8; 128];
        bytes[..7].copy_from_slice(b"NOTOURS");
        std::fs::write(&path, &bytes).unwrap();
        let err = MmapDb::open(&path).expect_err("must reject");
        assert!(matches!(err, MmapDbError::BadMagic), "got {err:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
