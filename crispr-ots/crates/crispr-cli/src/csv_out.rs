//! GuideScan2-compatible CSV writer for the `enumerate` subcommand.
//!
//! Matches the format that crisprware's `score_guides.py:217-237` reads
//! back:
//!
//! ```text
//! id,sequence,match_chrm,match_position,match_strand,match_distance,specificity
//! chr1:916690:+,CTGGCAGCATGAGCGATGGCNGG,chr11_74452305:74482685,5347,+,0,0.398429
//! ```
//!
//! One row per (guide, off-target hit) — including the on-target row
//! (`match_distance` = 0 at the guide's own coordinates) and any
//! multi-mapping copies. crisprware deduplicates by `id` and only keeps
//! the `specificity` column, so the per-hit shape is a free win for the
//! "report all off-targets up to N" use case the user asked for.

use std::io::{self, Write};

/// One CSV row.
#[derive(Debug, Clone)]
pub struct CsvRow<'a> {
    /// Caller-supplied guide id (echoed verbatim from the kmers CSV).
    pub id: &'a str,
    /// Guide sequence + PAM, concatenated as printed (e.g. with the
    /// literal "N" if the input PAM was NGG).
    pub sequence: &'a str,
    /// Off-target contig name.
    pub match_chrm: &'a str,
    /// 0-indexed off-target start on the contig.
    pub match_position: u32,
    /// `+` or `-`.
    pub match_strand: char,
    /// Mismatch count for this off-target row.
    pub match_distance: u8,
    /// Per-guide specificity, repeated on every row of the same guide id.
    /// 6 significant figures match GuideScan2's print convention; the
    /// downstream consumer (`pd.read_csv` in crisprware) doesn't care.
    pub specificity: f64,
}

pub const HEADER: &str = "id,sequence,match_chrm,match_position,match_strand,match_distance,specificity";

/// Write the CSV header plus all rows to `writer`.
///
/// # Errors
/// Surfaces any I/O error from the underlying writer.
pub fn write_rows<W: Write>(writer: &mut W, rows: &[CsvRow<'_>]) -> io::Result<()> {
    writeln!(writer, "{HEADER}")?;
    for row in rows {
        writeln!(
            writer,
            "{id},{sequence},{chrm},{pos},{strand},{dist},{spec:.6}",
            id = row.id,
            sequence = row.sequence,
            chrm = row.match_chrm,
            pos = row.match_position,
            strand = row.match_strand,
            dist = row.match_distance,
            spec = row.specificity,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_expected_header() {
        let mut buf = Vec::new();
        write_rows(&mut buf, &[]).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(
            text.trim_end(),
            "id,sequence,match_chrm,match_position,match_strand,match_distance,specificity"
        );
    }

    #[test]
    fn six_sigfig_specificity() {
        let mut buf = Vec::new();
        write_rows(
            &mut buf,
            &[CsvRow {
                id: "g1",
                sequence: "ACGTACGTACGTACGTACGTNGG",
                match_chrm: "chr1",
                match_position: 100,
                match_strand: '+',
                match_distance: 0,
                specificity: 0.398_428_5,
            }],
        )
        .unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data = text.lines().nth(1).unwrap();
        assert!(
            data.ends_with(",0.398429") || data.ends_with(",0.398428"),
            "expected 6-sigfig spec, got {data:?}"
        );
    }
}
