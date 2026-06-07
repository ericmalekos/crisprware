//! Reader for the per-guide CSV input format accepted by GuideScan2's
//! `enumerate --kmers-file` and produced by crisprware's `score_guides`.
//!
//! The file is a single-column TSV whose header is the literal string
//! `id,sequence,pam,chromosome,position,sense` (the comma in the column
//! name is what makes this a "kmers CSV" — pandas writes it with
//! `sep="\t"` but the single column already has comma-separated subfields).
//! Body rows are comma-separated values of that shape:
//!
//! ```text
//! id,sequence,pam,chromosome,position,sense
//! chr21:5241815:+,ATAACATTTCAGCATTTCCA,NGG,chr21,5241815,+
//! ```
//!
//! Fields:
//! - `id` — opaque identifier echoed back on every output row. Convention
//!   `<chrom>:<pos>:<strand>` but we treat it as a string.
//! - `sequence` — protospacer only (20 nt for SpCas9). No PAM.
//! - `pam` — PAM as listed (e.g. `NGG`). May contain IUPAC ambiguity
//!   codes; we treat any non-ACGT character as `A` because PAM bits are
//!   masked out of the off-target compare and not read for CFD scoring
//!   (CFD reads the *off-target's* PAM, not the guide's).
//! - `chromosome`, `position`, `sense` — on-target genomic coordinates.
//!   `position` is the 1-indexed start in GuideScan2's convention; we
//!   store it verbatim because the only consumer (the CSV writer) emits
//!   it back unchanged.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crispr_encoding::{Base, Enzyme, PamSide, Site};

/// One parsed entry from a kmers CSV. Owned, so the file handle can be
/// dropped right after parsing.
#[derive(Debug, Clone)]
pub struct KmerEntry {
    /// Caller-supplied identifier, echoed back on every output row.
    pub id: String,
    /// Protospacer bases as written (no PAM).
    pub sequence: String,
    /// PAM bases as written. May include IUPAC ambiguity codes.
    pub pam: String,
    /// On-target chromosome name, e.g. "chr21".
    pub chromosome: String,
    /// 1-indexed on-target start, in GuideScan2's convention.
    pub position: u64,
    /// On-target strand: `'+'` or `'-'`.
    pub sense: char,
    /// Bit-encoded full scan target (protospacer + PAM, ordered by the
    /// enzyme's `pam_side`). Non-ACGT bases in the PAM are encoded as
    /// `A` since they're masked out during off-target compare.
    pub site: Site,
}

/// Errors from parsing a kmers CSV.
#[derive(Debug)]
pub enum KmerCsvError {
    Io(std::io::Error),
    /// The first line didn't match the expected header.
    BadHeader {
        line: String,
    },
    /// A data row had the wrong number of comma-separated fields.
    BadRowArity {
        line_number: usize,
        expected: usize,
        actual: usize,
        line: String,
    },
    /// A row's `sequence` or `pam` field had a length that doesn't match
    /// the configured enzyme.
    BadSequenceLength {
        line_number: usize,
        kind: &'static str,
        expected: usize,
        actual: usize,
        line: String,
    },
    /// A row's `sequence` had a non-ACGT base.
    BadSequenceBase {
        line_number: usize,
        line: String,
    },
    /// A row's `position` field wasn't parseable as a `u64`.
    BadPosition {
        line_number: usize,
        raw: String,
    },
    /// A row's `sense` wasn't `+` or `-`.
    BadSense {
        line_number: usize,
        raw: String,
    },
    /// The file had no data rows.
    Empty,
}

impl std::fmt::Display for KmerCsvError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error reading kmers CSV: {e}"),
            Self::BadHeader { line } => write!(
                f,
                "kmers CSV header mismatch — expected \
                 `id,sequence,pam,chromosome,position,sense` (single column \
                 of comma-separated subfields), got: {line:?}"
            ),
            Self::BadRowArity {
                line_number,
                expected,
                actual,
                line,
            } => write!(
                f,
                "kmers CSV line {line_number}: expected {expected} comma-separated \
                 fields, got {actual}: {line:?}"
            ),
            Self::BadSequenceLength {
                line_number,
                kind,
                expected,
                actual,
                line,
            } => write!(
                f,
                "kmers CSV line {line_number}: {kind} length {actual} doesn't \
                 match enzyme's {kind} length {expected}: {line:?}"
            ),
            Self::BadSequenceBase { line_number, line } => write!(
                f,
                "kmers CSV line {line_number}: non-ACGT base in `sequence` \
                 field (PAM may contain IUPAC; protospacer may not): {line:?}"
            ),
            Self::BadPosition { line_number, raw } => write!(
                f,
                "kmers CSV line {line_number}: invalid `position` field \
                 (expected non-negative integer): {raw:?}"
            ),
            Self::BadSense { line_number, raw } => write!(
                f,
                "kmers CSV line {line_number}: invalid `sense` field \
                 (expected `+` or `-`): {raw:?}"
            ),
            Self::Empty => write!(f, "kmers CSV contained no data rows"),
        }
    }
}

impl std::error::Error for KmerCsvError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for KmerCsvError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

/// Expected single column header, written verbatim by crisprware's
/// `score_guides.py` and accepted by GuideScan2.
pub const HEADER: &str = "id,sequence,pam,chromosome,position,sense";

/// Parse the kmers CSV at `path` against the given `enzyme`.
///
/// # Errors
/// Returns [`KmerCsvError`] for I/O failures, malformed header/rows, or
/// an empty body.
pub fn read_from_path(path: &Path, enzyme: &Enzyme) -> Result<Vec<KmerEntry>, KmerCsvError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    read_from_reader(reader, enzyme)
}

/// Parse a kmers CSV from an arbitrary `BufRead`. Useful for tests that
/// want to feed in-memory bytes without touching the filesystem.
///
/// # Errors
/// Same as [`read_from_path`].
pub fn read_from_reader<R: BufRead>(
    reader: R,
    enzyme: &Enzyme,
) -> Result<Vec<KmerEntry>, KmerCsvError> {
    let mut lines = reader.lines();

    let header_line = match lines.next() {
        Some(line) => line?,
        None => return Err(KmerCsvError::Empty),
    };
    let header_trimmed = header_line.trim_end_matches(['\r', '\n']).trim();
    if header_trimmed != HEADER {
        return Err(KmerCsvError::BadHeader { line: header_line });
    }

    let proto_len = usize::from(enzyme.protospacer_len);
    let pam_len = usize::from(enzyme.pam_len);
    let total_len = proto_len + pam_len;

    let mut entries: Vec<KmerEntry> = Vec::new();
    for (idx, line_result) in lines.enumerate() {
        let raw = line_result?;
        let line = raw.trim_end_matches(['\r', '\n']);
        if line.trim().is_empty() {
            continue;
        }
        let line_number = idx + 2; // header is line 1; first body row is line 2
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != 6 {
            return Err(KmerCsvError::BadRowArity {
                line_number,
                expected: 6,
                actual: fields.len(),
                line: raw,
            });
        }

        let id = fields[0].to_string();
        let sequence = fields[1].to_string();
        let pam = fields[2].to_string();
        let chromosome = fields[3].to_string();
        let position: u64 = fields[4].parse().map_err(|_| KmerCsvError::BadPosition {
            line_number,
            raw: fields[4].to_string(),
        })?;
        let sense = match fields[5].trim() {
            "+" => '+',
            "-" => '-',
            other => {
                return Err(KmerCsvError::BadSense {
                    line_number,
                    raw: other.to_string(),
                });
            }
        };

        if sequence.len() != proto_len {
            return Err(KmerCsvError::BadSequenceLength {
                line_number,
                kind: "protospacer",
                expected: proto_len,
                actual: sequence.len(),
                line: raw,
            });
        }
        if pam.len() != pam_len {
            return Err(KmerCsvError::BadSequenceLength {
                line_number,
                kind: "PAM",
                expected: pam_len,
                actual: pam.len(),
                line: raw,
            });
        }

        let site = encode_site(&sequence, &pam, enzyme.pam_side, total_len).ok_or_else(|| {
            KmerCsvError::BadSequenceBase {
                line_number,
                line: raw.clone(),
            }
        })?;

        entries.push(KmerEntry {
            id,
            sequence,
            pam,
            chromosome,
            position,
            sense,
            site,
        });
    }

    if entries.is_empty() {
        return Err(KmerCsvError::Empty);
    }
    Ok(entries)
}

/// Encode protospacer + PAM into a [`Site`]. The two are concatenated in
/// the order dictated by `pam_side` (3' PAM → protospacer first; 5' PAM →
/// PAM first). Non-ACGT bases in the protospacer cause a `None` return;
/// non-ACGT bases in the PAM are silently encoded as `A` (PAM bits are
/// masked out of the off-target compare and not read for CFD scoring).
fn encode_site(protospacer: &str, pam: &str, side: PamSide, total_len: usize) -> Option<Site> {
    let mut buf: Vec<Base> = Vec::with_capacity(total_len);
    let push_protospacer = |buf: &mut Vec<Base>, p: &str| -> Option<()> {
        for c in p.bytes() {
            buf.push(base_from_ascii(c)?);
        }
        Some(())
    };
    let push_pam = |buf: &mut Vec<Base>, p: &str| {
        for c in p.bytes() {
            buf.push(base_from_ascii(c).unwrap_or(Base::A));
        }
    };
    match side {
        PamSide::ThreePrime => {
            push_protospacer(&mut buf, protospacer)?;
            push_pam(&mut buf, pam);
        }
        PamSide::FivePrime => {
            push_pam(&mut buf, pam);
            push_protospacer(&mut buf, protospacer)?;
        }
    }
    debug_assert_eq!(buf.len(), total_len);
    Some(Site::encode(buf))
}

fn base_from_ascii(c: u8) -> Option<Base> {
    match c.to_ascii_uppercase() {
        b'A' => Some(Base::A),
        b'C' => Some(Base::C),
        b'G' => Some(Base::G),
        b'T' => Some(Base::T),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_csv(body: &str) -> std::io::Cursor<Vec<u8>> {
        std::io::Cursor::new(body.as_bytes().to_vec())
    }

    #[test]
    fn parses_two_row_fixture() {
        let body = "\
id,sequence,pam,chromosome,position,sense
chr21:5241815:+,ATAACATTTCAGCATTTCCA,NGG,chr21,5241815,+
chr21:46681562:-,GAATGTGCTGATAGCGGACG,NGG,chr21,46681562,-
";
        let enzyme = Enzyme::spcas9_ngg();
        let parsed =
            read_from_reader(BufReader::new(write_csv(body)), &enzyme).expect("parse succeeds");
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].id, "chr21:5241815:+");
        assert_eq!(parsed[0].sequence, "ATAACATTTCAGCATTTCCA");
        assert_eq!(parsed[0].pam, "NGG");
        assert_eq!(parsed[0].chromosome, "chr21");
        assert_eq!(parsed[0].position, 5_241_815);
        assert_eq!(parsed[0].sense, '+');
        assert_eq!(parsed[1].sense, '-');
    }

    #[test]
    fn rejects_bad_header() {
        let body = "wrong header\nrow,row,row,row,row,row\n";
        let enzyme = Enzyme::spcas9_ngg();
        match read_from_reader(BufReader::new(write_csv(body)), &enzyme) {
            Err(KmerCsvError::BadHeader { .. }) => {}
            other => panic!("expected BadHeader, got {other:?}"),
        }
    }

    #[test]
    fn rejects_wrong_protospacer_length() {
        let body = "\
id,sequence,pam,chromosome,position,sense
short,ACGT,NGG,chr1,100,+
";
        let enzyme = Enzyme::spcas9_ngg();
        match read_from_reader(BufReader::new(write_csv(body)), &enzyme) {
            Err(KmerCsvError::BadSequenceLength { kind, .. }) => assert_eq!(kind, "protospacer"),
            other => panic!("expected BadSequenceLength, got {other:?}"),
        }
    }

    #[test]
    fn rejects_non_acgt_in_protospacer() {
        let body = "\
id,sequence,pam,chromosome,position,sense
g,NTAACATTTCAGCATTTCCA,NGG,chr1,1,+
";
        let enzyme = Enzyme::spcas9_ngg();
        match read_from_reader(BufReader::new(write_csv(body)), &enzyme) {
            Err(KmerCsvError::BadSequenceBase { .. }) => {}
            other => panic!("expected BadSequenceBase, got {other:?}"),
        }
    }

    #[test]
    fn accepts_iupac_in_pam_field() {
        // `NGG` is the typical value crisprware writes; the N should be
        // accepted and silently encoded as A.
        let body = "\
id,sequence,pam,chromosome,position,sense
g,ATAACATTTCAGCATTTCCA,NGG,chr1,1,+
";
        let enzyme = Enzyme::spcas9_ngg();
        let parsed =
            read_from_reader(BufReader::new(write_csv(body)), &enzyme).expect("parse succeeds");
        assert_eq!(parsed.len(), 1);
        let decoded = parsed[0]
            .site
            .decode_ascii(usize::from(enzyme.total_scan_len()));
        // Protospacer is intact; PAM `N` was substituted by `A`.
        assert_eq!(decoded, "ATAACATTTCAGCATTTCCAAGG");
    }

    #[test]
    fn empty_body_is_an_error() {
        let body = "id,sequence,pam,chromosome,position,sense\n";
        let enzyme = Enzyme::spcas9_ngg();
        assert!(matches!(
            read_from_reader(BufReader::new(write_csv(body)), &enzyme),
            Err(KmerCsvError::Empty)
        ));
    }
}
