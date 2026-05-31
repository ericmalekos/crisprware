//! Minimal streaming FASTA reader.
//!
//! Handles plain `.fa`/`.fasta` and gzip-compressed `.fa.gz` transparently
//! based on file extension. Sequence bytes are upper-cased; whitespace
//! within sequence lines is stripped. Non-ACGT bases (Ns, IUPAC ambiguity)
//! are preserved as-is — the site finder filters them out on a per-window
//! basis.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

/// One contig: name + uppercase ASCII sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub sequence: Vec<u8>,
}

/// Read every contig from a file. `.gz` extension triggers gzip decoding.
pub fn read_fasta(path: impl AsRef<Path>) -> io::Result<Vec<Contig>> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let reader: Box<dyn Read> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    read_fasta_from(BufReader::new(reader))
}

/// Read every contig from any `BufRead` source.
pub fn read_fasta_from<R: BufRead>(mut reader: R) -> io::Result<Vec<Contig>> {
    let mut contigs = Vec::new();
    let mut current: Option<Contig> = None;
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.is_empty() {
            continue;
        }
        if let Some(stripped) = trimmed.strip_prefix('>') {
            if let Some(c) = current.take() {
                contigs.push(c);
            }
            // FASTA header convention: contig name is the first whitespace-
            // delimited token after the `>`. Anything after that is treated
            // as a description and dropped.
            let name = stripped.split_whitespace().next().unwrap_or("").to_string();
            current = Some(Contig {
                name,
                sequence: Vec::new(),
            });
        } else if let Some(c) = current.as_mut() {
            for b in trimmed.bytes() {
                if !b.is_ascii_whitespace() {
                    c.sequence.push(b.to_ascii_uppercase());
                }
            }
        }
        // Otherwise: sequence line before any `>` header — silently ignored.
    }
    if let Some(c) = current {
        contigs.push(c);
    }
    Ok(contigs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    #[test]
    fn parse_two_contigs() {
        let data = b">chr1 some description\nACGT\nACGT\n>chr2\nTTTT\n";
        let contigs = read_fasta_from(&data[..]).expect("parses");
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[0].name, "chr1");
        assert_eq!(contigs[0].sequence, b"ACGTACGT");
        assert_eq!(contigs[1].name, "chr2");
        assert_eq!(contigs[1].sequence, b"TTTT");
    }

    #[test]
    fn lowercase_uppercased_and_blank_lines_skipped() {
        let data = b">chr1\nacgt\n\nacgt\n";
        let contigs = read_fasta_from(&data[..]).expect("parses");
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].sequence, b"ACGTACGT");
    }

    #[test]
    fn ns_are_preserved() {
        let data = b">chr1\nACGTNNNACGT\n";
        let contigs = read_fasta_from(&data[..]).expect("parses");
        assert_eq!(contigs[0].sequence, b"ACGTNNNACGT");
    }

    #[test]
    fn crlf_line_endings() {
        let data = b">chr1\r\nACGT\r\nACGT\r\n";
        let contigs = read_fasta_from(&data[..]).expect("parses");
        assert_eq!(contigs[0].sequence, b"ACGTACGT");
    }

    #[test]
    fn empty_input_yields_empty_list() {
        let data: &[u8] = b"";
        let contigs = read_fasta_from(data).expect("parses");
        assert!(contigs.is_empty());
    }

    #[test]
    fn gzipped_roundtrip() {
        // Build a temp file with gzipped FASTA content, then read it back
        // through the dispatch path.
        let dir = tempdir_in_target();
        let path = dir.join("test.fa.gz");
        let plain = b">chr1\nACGT\nACGT\n";
        let file = File::create(&path).expect("creates temp file");
        let mut enc = GzEncoder::new(file, Compression::default());
        enc.write_all(plain).expect("gzip writes");
        enc.finish().expect("gzip finishes");

        let contigs = read_fasta(&path).expect("parses");
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].sequence, b"ACGTACGT");

        let _ = std::fs::remove_file(&path);
    }

    /// Allocate a unique temp directory under `target/` so tests on shared
    /// filesystems don't collide. Created lazily; not cleaned up
    /// automatically (the workspace's `.gitignore` ignores `target/`).
    fn tempdir_in_target() -> std::path::PathBuf {
        let base = std::env::temp_dir();
        let name = format!(
            "crispr-ots-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_nanos())
        );
        let path = base.join(name);
        std::fs::create_dir_all(&path).expect("creates temp dir");
        path
    }
}
