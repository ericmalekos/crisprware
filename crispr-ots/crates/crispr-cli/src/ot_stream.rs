//! Streaming writers for the scalable scored-output modes (`--output-mode`).
//!
//! Two outputs, both derived from one streaming pass over the GPU hit batches
//! (see `discover::run_discover_streaming`):
//!
//! * **Mode 1 (aggregated):** one row per guide — `id, specificity[, max_cfd]`
//!   ([`AggWriter`]). Computed over *all* off-targets (no CFD floor).
//! * **Mode 2 (per-off-target):** one record per (guide, off-target) —
//!   `guide_id, contig, start, strand, CFD` — written either as a compact
//!   12-byte binary record ([`OtBinaryWriter`]) or TSV ([`OtTsvWriter`]),
//!   floored at `--cfd-threshold` for disk control. A [`GuidesSidecarWriter`]
//!   maps `guide_id → (id, on-target locus, sequence, on_count)` so binary
//!   records stay tiny but joinable.
//!
//! Binary layout — 64-byte header then `n` × 12-byte records:
//! ```text
//! header: magic "CROTOT2\0" (8) | version u16 | record_bytes u16 | max_mm u8 |
//!         (pad 3) | cfd_threshold f64 | n_guides u64 | (pad to 64)
//! record: guide_id u32 | contig u16 | offset(bits0-30)+strand(bit31) u32 | cfd_q u16
//! ```
//! CFD is quantized to `u16` as `round(cfd * 65535)` (precision ~1.5e-5, far
//! below biological significance); strand is folded into bit 31 of the offset
//! (mirroring `PackedEntry`'s own convention).

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Magic bytes opening a Mode-2 binary off-target file.
pub const OT_MAGIC: [u8; 8] = *b"CROTOT2\0";
/// Bytes per binary off-target record.
pub const OT_RECORD_BYTES: usize = 12;
/// Fixed binary header size.
pub const OT_HEADER_BYTES: usize = 64;
/// Current binary format version.
pub const OT_VERSION: u16 = 1;
const CFD_SCALE: f64 = 65535.0;

/// One per-off-target record (12 bytes LE).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OtRecord {
    /// Absolute guide index (joins to the sidecar / Mode-1 row order).
    pub guide_id: u32,
    /// Off-target contig id (`Position::contig_id`, ≤ u16 by index invariant).
    pub contig: u16,
    /// Off-target 0-based genomic offset (low 31 bits used).
    pub offset: u32,
    /// `true` if the off-target is on the reverse strand.
    pub strand_reverse: bool,
    /// CFD quantized to `u16` (`round(cfd * 65535)`).
    pub cfd_q: u16,
}

impl OtRecord {
    /// Quantize a CFD score in `[0, 1]` to `u16`.
    #[must_use]
    pub fn quantize_cfd(cfd: f64) -> u16 {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let q = (cfd.clamp(0.0, 1.0) * CFD_SCALE).round() as u16;
        q
    }

    /// Recover the (lossy) CFD score from a quantized value.
    #[must_use]
    pub fn dequantize_cfd(q: u16) -> f64 {
        f64::from(q) / CFD_SCALE
    }

    /// This record's dequantized CFD.
    #[must_use]
    pub fn cfd(self) -> f64 {
        Self::dequantize_cfd(self.cfd_q)
    }

    #[must_use]
    fn to_bytes(self) -> [u8; OT_RECORD_BYTES] {
        let off_strand = (self.offset & 0x7FFF_FFFF) | (u32::from(self.strand_reverse) << 31);
        let mut b = [0u8; OT_RECORD_BYTES];
        b[0..4].copy_from_slice(&self.guide_id.to_le_bytes());
        b[4..6].copy_from_slice(&self.contig.to_le_bytes());
        b[6..10].copy_from_slice(&off_strand.to_le_bytes());
        b[10..12].copy_from_slice(&self.cfd_q.to_le_bytes());
        b
    }

    /// Decode a 12-byte little-endian record.
    #[must_use]
    pub fn from_bytes(b: &[u8; OT_RECORD_BYTES]) -> Self {
        let guide_id = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        let contig = u16::from_le_bytes([b[4], b[5]]);
        let off_strand = u32::from_le_bytes([b[6], b[7], b[8], b[9]]);
        let cfd_q = u16::from_le_bytes([b[10], b[11]]);
        Self {
            guide_id,
            contig,
            offset: off_strand & 0x7FFF_FFFF,
            strand_reverse: (off_strand >> 31) == 1,
            cfd_q,
        }
    }
}

/// Self-describing header for the Mode-2 binary file.
#[derive(Debug, Clone, Copy)]
pub struct OtHeader {
    /// Max mismatches the scan used.
    pub max_mm: u8,
    /// CFD floor applied to the written records (`0.0` = none).
    pub cfd_threshold: f64,
    /// Number of guides scored.
    pub n_guides: u64,
}

impl OtHeader {
    #[must_use]
    fn to_bytes(self) -> [u8; OT_HEADER_BYTES] {
        let mut b = [0u8; OT_HEADER_BYTES];
        b[0..8].copy_from_slice(&OT_MAGIC);
        b[8..10].copy_from_slice(&OT_VERSION.to_le_bytes());
        #[allow(clippy::cast_possible_truncation)]
        let rb = OT_RECORD_BYTES as u16;
        b[10..12].copy_from_slice(&rb.to_le_bytes());
        b[12] = self.max_mm;
        b[16..24].copy_from_slice(&self.cfd_threshold.to_le_bytes());
        b[24..32].copy_from_slice(&self.n_guides.to_le_bytes());
        b
    }

    /// Read and validate a header from `r`.
    ///
    /// # Errors
    /// Returns an error on I/O failure or a bad magic.
    pub fn read_from<R: Read>(r: &mut R) -> io::Result<Self> {
        let mut b = [0u8; OT_HEADER_BYTES];
        r.read_exact(&mut b)?;
        if b[0..8] != OT_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "not a CROTOT2 off-target file (bad magic)",
            ));
        }
        Ok(Self {
            max_mm: b[12],
            cfd_threshold: f64::from_le_bytes(b[16..24].try_into().expect("8 bytes")),
            n_guides: u64::from_le_bytes(b[24..32].try_into().expect("8 bytes")),
        })
    }
}

/// Mode-2 binary off-target writer.
#[derive(Debug)]
pub struct OtBinaryWriter {
    w: BufWriter<File>,
}

impl OtBinaryWriter {
    /// Create the file and write the header.
    ///
    /// # Errors
    /// I/O failure creating or writing the file.
    pub fn create(path: &Path, header: OtHeader) -> io::Result<Self> {
        let mut w = BufWriter::with_capacity(1 << 20, File::create(path)?);
        w.write_all(&header.to_bytes())?;
        Ok(Self { w })
    }

    /// Append one record.
    ///
    /// # Errors
    /// I/O failure.
    pub fn write_record(&mut self, rec: OtRecord) -> io::Result<()> {
        self.w.write_all(&rec.to_bytes())
    }

    /// Append a pre-encoded byte chunk (e.g. concatenated thread-local buffers).
    ///
    /// # Errors
    /// I/O failure.
    pub fn write_chunk(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.w.write_all(bytes)
    }

    /// Flush and close.
    ///
    /// # Errors
    /// I/O failure.
    pub fn finish(mut self) -> io::Result<()> {
        self.w.flush()
    }
}

/// Encode a record into a caller-owned buffer (for parallel batch encoding).
pub fn push_record_bytes(buf: &mut Vec<u8>, rec: OtRecord) {
    buf.extend_from_slice(&rec.to_bytes());
}

/// Mode-2 TSV off-target writer (`guide_id\tchrom\tstart\tstrand\tcfd`).
#[derive(Debug)]
pub struct OtTsvWriter {
    w: BufWriter<File>,
}

impl OtTsvWriter {
    /// Create the file and write the header line.
    ///
    /// # Errors
    /// I/O failure.
    pub fn create(path: &Path) -> io::Result<Self> {
        let mut w = BufWriter::with_capacity(1 << 20, File::create(path)?);
        writeln!(w, "guide_id\tchrom\tstart\tstrand\tcfd")?;
        Ok(Self { w })
    }

    /// Append one row.
    ///
    /// # Errors
    /// I/O failure.
    pub fn write_row(
        &mut self,
        guide_id: u32,
        chrom: &str,
        start: u32,
        strand_reverse: bool,
        cfd: f64,
    ) -> io::Result<()> {
        let strand = if strand_reverse { '-' } else { '+' };
        writeln!(self.w, "{guide_id}\t{chrom}\t{start}\t{strand}\t{cfd:.6}")
    }

    /// Flush and close.
    ///
    /// # Errors
    /// I/O failure.
    pub fn finish(mut self) -> io::Result<()> {
        self.w.flush()
    }
}

/// Mode-1 aggregated per-guide writer.
#[derive(Debug)]
pub struct AggWriter {
    w: BufWriter<File>,
}

impl AggWriter {
    /// Create the file and write the header. `cas12a` selects the column set.
    ///
    /// # Errors
    /// I/O failure.
    pub fn create(path: &Path, cas12a: bool) -> io::Result<Self> {
        let mut w = BufWriter::new(File::create(path)?);
        if cas12a {
            writeln!(
                w,
                "id,specificity_tttn,specificity_tttv,max_cfd,off_target_count,saturated"
            )?;
        } else {
            writeln!(w, "id,specificity,max_cfd,off_target_count,saturated")?;
        }
        Ok(Self { w })
    }

    /// SpCas9 row. `saturated` marks a guide whose off-target cap was hit, so
    /// `specificity`/`off_target_count` are partial (specificity is then an
    /// upper bound — the guide is effectively non-specific).
    ///
    /// # Errors
    /// I/O failure.
    pub fn write_cas9(
        &mut self,
        id: &str,
        specificity: f64,
        max_cfd: f64,
        off_target_count: u32,
        saturated: bool,
    ) -> io::Result<()> {
        writeln!(
            self.w,
            "{id},{specificity:.6},{max_cfd:.6},{off_target_count},{}",
            u8::from(saturated)
        )
    }

    /// Cas12a row. See [`write_cas9`](Self::write_cas9) for `saturated`.
    ///
    /// # Errors
    /// I/O failure.
    pub fn write_cas12a(
        &mut self,
        id: &str,
        tttn: f64,
        tttv: f64,
        max_cfd: f64,
        off_target_count: u32,
        saturated: bool,
    ) -> io::Result<()> {
        writeln!(
            self.w,
            "{id},{tttn:.6},{tttv:.6},{max_cfd:.6},{off_target_count},{}",
            u8::from(saturated)
        )
    }

    /// Flush and close.
    ///
    /// # Errors
    /// I/O failure.
    pub fn finish(mut self) -> io::Result<()> {
        self.w.flush()
    }
}

/// Sidecar mapping `guide_id` → on-target metadata for Mode-2 joins.
#[derive(Debug)]
pub struct GuidesSidecarWriter {
    w: BufWriter<File>,
}

impl GuidesSidecarWriter {
    /// Create the sidecar and write its header.
    ///
    /// # Errors
    /// I/O failure.
    pub fn create(path: &Path) -> io::Result<Self> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "guide_id\tid\tchrom\tposition\tstrand\tsequence\ton_count")?;
        Ok(Self { w })
    }

    /// One guide row.
    ///
    /// # Errors
    /// I/O failure.
    #[allow(clippy::too_many_arguments)]
    pub fn write_row(
        &mut self,
        guide_id: u32,
        id: &str,
        chrom: &str,
        position: u32,
        strand_reverse: bool,
        sequence: &str,
        on_count: u32,
    ) -> io::Result<()> {
        let strand = if strand_reverse { '-' } else { '+' };
        writeln!(
            self.w,
            "{guide_id}\t{id}\t{chrom}\t{position}\t{strand}\t{sequence}\t{on_count}"
        )
    }

    /// Flush and close.
    ///
    /// # Errors
    /// I/O failure.
    pub fn finish(mut self) -> io::Result<()> {
        self.w.flush()
    }
}

/// Per-guide aggregate reproduced from a Mode-2 binary file.
#[derive(Debug, Clone, Copy)]
pub struct AggFromOt {
    /// FlashFry specificity `1 / (1 + Σ cfd)` over the file's records.
    pub specificity: f64,
    /// Largest CFD seen for the guide.
    pub max_cfd: f64,
    /// Number of off-target records for the guide.
    pub n_ot: u64,
}

/// Stream a Mode-2 binary file and reproduce each guide's FlashFry
/// specificity (`1 / (1 + Σ cfd)`). The validation tool for "Mode 2 → Mode 1";
/// with `--cfd-threshold 0` (no floor) this equals the Mode-1 output exactly
/// (on the FlashFry convention; off-targets only — on-targets aren't written).
///
/// # Errors
/// I/O failure or a malformed file.
pub fn aggregate_binary_flashfry(path: &Path) -> io::Result<HashMap<u32, AggFromOt>> {
    let mut r = BufReader::with_capacity(1 << 20, File::open(path)?);
    let _hdr = OtHeader::read_from(&mut r)?;
    let mut acc: HashMap<u32, (f64, f64, u64)> = HashMap::new();
    let mut buf = [0u8; OT_RECORD_BYTES];
    loop {
        match r.read_exact(&mut buf) {
            Ok(()) => {
                let rec = OtRecord::from_bytes(&buf);
                let cfd = rec.cfd();
                let e = acc.entry(rec.guide_id).or_insert((0.0, 0.0, 0));
                e.0 += cfd;
                if cfd > e.1 {
                    e.1 = cfd;
                }
                e.2 += 1;
            }
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
    }
    Ok(acc
        .into_iter()
        .map(|(g, (sum, max, n))| {
            (
                g,
                AggFromOt {
                    specificity: 1.0 / (1.0 + sum),
                    max_cfd: max,
                    n_ot: n,
                },
            )
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn record_round_trips() {
        let rec = OtRecord {
            guide_id: 123_456,
            contig: 21,
            offset: 0x4321_0FED & 0x7FFF_FFFF,
            strand_reverse: true,
            cfd_q: OtRecord::quantize_cfd(0.37),
        };
        let bytes = rec.to_bytes();
        assert_eq!(bytes.len(), OT_RECORD_BYTES);
        let back = OtRecord::from_bytes(&bytes);
        assert_eq!(rec, back);
        assert!((back.cfd() - 0.37).abs() < 1e-4);
    }

    #[test]
    fn cfd_quantization_bounds() {
        assert_eq!(OtRecord::quantize_cfd(0.0), 0);
        assert_eq!(OtRecord::quantize_cfd(1.0), 65535);
        assert_eq!(OtRecord::quantize_cfd(2.0), 65535); // clamped
        assert_eq!(OtRecord::quantize_cfd(-1.0), 0); // clamped
    }

    #[test]
    fn header_round_trips() {
        let hdr = OtHeader {
            max_mm: 4,
            cfd_threshold: 0.023,
            n_guides: 5_295_011,
        };
        let bytes = hdr.to_bytes();
        assert_eq!(bytes.len(), OT_HEADER_BYTES);
        let mut cur = Cursor::new(bytes.to_vec());
        let back = OtHeader::read_from(&mut cur).expect("reads back");
        assert_eq!(back.max_mm, 4);
        assert!((back.cfd_threshold - 0.023).abs() < 1e-12);
        assert_eq!(back.n_guides, 5_295_011);
    }

    #[test]
    fn aggregator_reproduces_flashfry_specificity() {
        // Two guides; guide 0 has two off-targets (cfd 0.5, 0.25), guide 1 one (0.1).
        let dir = std::env::temp_dir();
        let path = dir.join(format!("crispr_ot_test_{}.bin", std::process::id()));
        let mut w = OtBinaryWriter::create(
            &path,
            OtHeader {
                max_mm: 4,
                cfd_threshold: 0.0,
                n_guides: 2,
            },
        )
        .expect("create");
        for (g, cfd) in [(0u32, 0.5_f64), (0, 0.25), (1, 0.1)] {
            w.write_record(OtRecord {
                guide_id: g,
                contig: 0,
                offset: 100,
                strand_reverse: false,
                cfd_q: OtRecord::quantize_cfd(cfd),
            })
            .expect("write");
        }
        w.finish().expect("finish");

        let agg = aggregate_binary_flashfry(&path).expect("aggregate");
        // guide 0: 1/(1 + 0.5 + 0.25) = 1/1.75
        assert!((agg[&0].specificity - 1.0 / 1.75).abs() < 1e-3);
        assert_eq!(agg[&0].n_ot, 2);
        assert!((agg[&0].max_cfd - 0.5).abs() < 1e-3);
        // guide 1: 1/(1 + 0.1)
        assert!((agg[&1].specificity - 1.0 / 1.1).abs() < 1e-3);
        std::fs::remove_file(&path).ok();
    }
}
