//! TSV writer for discover output.
//!
//! Column order is FlashFry-inspired but pared down to the fields the day-1
//! pipeline can produce. The optional `cfd_max` / `cfd_specificity` columns
//! appear only when `--score cfd` is requested, so a downstream consumer
//! can skip them if absent.

use std::io::{self, Write};

use crispr_db::Strand;

/// One discover result: a guide + its aggregated off-targets + (optionally)
/// scoring outputs.
#[derive(Debug, Clone)]
pub struct DiscoverRow {
    /// Contig name from the query FASTA where the guide was found.
    pub contig: String,
    /// 0-indexed start of the guide on the query contig (forward strand).
    pub start: u32,
    /// Exclusive end of the guide on the query contig.
    pub stop: u32,
    /// 23-bp (or enzyme-appropriate) protospacer+PAM in canonical order.
    pub target: String,
    /// Forward or reverse.
    pub strand: Strand,
    /// Number of *distinct* off-target sequences with mismatch count > 0.
    pub ot_count: usize,
    /// Comma-separated `SEQ_count_mismatches` entries, FlashFry-compatible.
    pub ot_sequences: String,
    /// Optional CFD outputs. `None` when `--score cfd` was not requested.
    pub cfd: Option<CfdColumns>,
}

#[derive(Debug, Clone, Copy)]
pub struct CfdColumns {
    pub cfd_max: f64,
    pub cfd_specificity: f64,
}

/// Write rows to a writer. The header includes CFD columns iff *any* row
/// has CFD data — the assumption is that all rows in a single run share
/// the same scoring configuration.
///
/// # Errors
/// Surfaces any I/O error from the underlying writer.
pub fn write_rows<W: Write>(writer: &mut W, rows: &[DiscoverRow]) -> io::Result<()> {
    let include_cfd = rows.iter().any(|r| r.cfd.is_some());
    if include_cfd {
        writeln!(
            writer,
            "contig\tstart\tstop\ttarget\torientation\totCount\totSequences\tcfd_max\tcfd_specificity"
        )?;
    } else {
        writeln!(
            writer,
            "contig\tstart\tstop\ttarget\torientation\totCount\totSequences"
        )?;
    }
    for row in rows {
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.contig,
            row.start,
            row.stop,
            row.target,
            row.strand.as_char(),
            row.ot_count,
            row.ot_sequences,
        )?;
        if include_cfd {
            if let Some(cfd) = row.cfd {
                write!(writer, "\t{}\t{}", cfd.cfd_max, cfd.cfd_specificity)?;
            } else {
                write!(writer, "\tNA\tNA")?;
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(target: &str, ot_count: usize) -> DiscoverRow {
        DiscoverRow {
            contig: "chr1".into(),
            start: 100,
            stop: 123,
            target: target.into(),
            strand: Strand::Forward,
            ot_count,
            ot_sequences: String::new(),
            cfd: None,
        }
    }

    #[test]
    fn no_cfd_header() {
        let rows = vec![row("AAAA...", 0)];
        let mut buf = Vec::new();
        write_rows(&mut buf, &rows).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let first_line = text.lines().next().unwrap();
        assert!(!first_line.contains("cfd"));
        assert!(first_line.contains("otCount"));
    }

    #[test]
    fn cfd_header_when_any_row_scored() {
        let mut row = row("AAAA...", 0);
        row.cfd = Some(CfdColumns {
            cfd_max: 0.1,
            cfd_specificity: 0.9,
        });
        let mut buf = Vec::new();
        write_rows(&mut buf, &[row]).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.lines().next().unwrap().contains("cfd_max"));
        assert!(text.lines().nth(1).unwrap().contains("0.1"));
    }

    #[test]
    fn strand_char_serialised() {
        let mut row = row("AAAA", 0);
        row.strand = Strand::Reverse;
        let mut buf = Vec::new();
        write_rows(&mut buf, &[row]).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // Data row has the strand char in column 5.
        let fields: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(fields[4], "-");
    }
}
