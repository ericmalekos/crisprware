//! End-to-end equivalence test against FlashFry on the chr22 quickstart
//! fixture. Auto-skips if the fixture is not present, so the test is safe
//! to leave on by default — set up `/tmp/ffry-quickstart` (see header
//! comment in `tests/chr22_equivalence_setup.md`-style instructions
//! emitted on skip) to actually exercise it.
//!
//! What the test validates:
//! 1. Our `SiteFinder` discovers the same set of guides in the EMX1 input
//!    FASTA as FlashFry's `discover` subcommand.
//! 2. For each guide, the count of *distinct* off-target sequences found
//!    in chr22 at up to 4 mismatches matches FlashFry's `otCount` column.
//! 3. Bonus: a per-guide off-target sequence set agrees one-for-one.

use std::collections::HashMap;
use std::collections::HashSet;
use std::path::Path;

use crispr_db::{read_fasta, BinTable, Position, SiteFinder, Strand};
use crispr_encoding::{Enzyme, Site};
use crispr_scan::{BinScanner, Guide, Scanner};

const CHR22_PATH: &str = "/tmp/ffry-quickstart/chr22.fa.gz";
const EMX1_PATH: &str = "/tmp/ffry-quickstart/EMX1_GAGTCCGAGCAGAAGAAGAAGGG.fasta";
const FLASHFRY_OUTPUT_PATH: &str = "/tmp/ffry-quickstart/EMX1.output";
/// FlashFry's `discover` default maximum mismatch count.
const MAX_MISMATCHES: u8 = 4;

/// One row of FlashFry's `EMX1.output`. We only parse the fields we need.
#[derive(Debug)]
struct FlashFryRow {
    /// 23-bp guide (protospacer + PAM) — column 4 ("target").
    target: String,
    /// Forward or reverse — column 7 ("orientation"), FWD or RVS.
    orientation: String,
    /// Total off-target count — column 8 ("otCount").
    ot_count: usize,
    /// Distinct off-target sequences with their per-OT mismatch counts.
    off_target_sequences: Vec<(String, u8)>,
}

fn parse_flashfry_output(path: &str) -> Vec<FlashFryRow> {
    let text = std::fs::read_to_string(path).expect("read FlashFry output");
    let mut lines = text.lines();
    let _header = lines.next();
    lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert!(
                fields.len() >= 9,
                "expected at least 9 tab-delimited fields, got {}: {line}",
                fields.len()
            );
            let target = fields[3].to_string();
            let orientation = fields[6].to_string();
            let ot_count: usize = fields[7].parse().expect("otCount is integer");
            // offTargets is a comma-separated list of "SEQ_count_mismatches"
            // entries. Parse the (SEQ, mismatches) pairs.
            let off_target_sequences: Vec<(String, u8)> = if fields[8].is_empty() {
                Vec::new()
            } else {
                fields[8]
                    .split(',')
                    .map(|entry| {
                        let parts: Vec<&str> = entry.rsplitn(3, '_').collect();
                        // rsplitn gives [mismatches, count, SEQ] in reverse order.
                        assert_eq!(parts.len(), 3, "malformed OT entry: {entry}");
                        let mismatches: u8 = parts[0].parse().expect("mm is integer");
                        let seq = parts[2].to_string();
                        (seq, mismatches)
                    })
                    .collect()
            };
            FlashFryRow {
                target,
                orientation,
                ot_count,
                off_target_sequences,
            }
        })
        .collect()
}

fn fixture_present() -> bool {
    Path::new(CHR22_PATH).exists()
        && Path::new(EMX1_PATH).exists()
        && Path::new(FLASHFRY_OUTPUT_PATH).exists()
}

/// Build a `BinTable` from a FASTA path using the given enzyme. Returns
/// the table plus the contig-id-to-name map (useful for reporting).
fn build_table(path: &str, enzyme: Enzyme) -> (BinTable, Vec<String>) {
    let contigs = read_fasta(path).expect("read FASTA");
    let names: Vec<String> = contigs.iter().map(|c| c.name.clone()).collect();
    let finder = SiteFinder::new(enzyme.clone());
    let mut table = BinTable::for_enzyme(enzyme);
    for (i, contig) in contigs.iter().enumerate() {
        finder.scan(
            &contig.sequence,
            u32::try_from(i).expect("contig fits in u32"),
            &mut table,
        );
    }
    (table, names)
}

fn discover_guides(path: &str, finder: &SiteFinder) -> Vec<(Site, Position)> {
    let contigs = read_fasta(path).expect("read EMX1 FASTA");
    let mut sites: Vec<(Site, Position)> = Vec::new();
    for (i, contig) in contigs.iter().enumerate() {
        finder.scan(
            &contig.sequence,
            u32::try_from(i).expect("contig fits in u32"),
            &mut sites,
        );
    }
    sites
}

#[test]
#[allow(clippy::too_many_lines)] // integration test, must run end-to-end
fn chr22_equivalence_against_flashfry() {
    if !fixture_present() {
        eprintln!(
            "Skipping chr22 equivalence test: fixture missing.\n\
             To set up:\n  \
             cd /tmp/ffry-quickstart\n  \
             tar xzf /path/to/FlashFry/test_data/quickstart_data.tar.gz\n  \
             wget https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar\n  \
             mkdir tmp\n  \
             java -Xmx4g -jar FlashFry-assembly-1.15.jar index \\\n    \
                 --tmpLocation ./tmp --database chr22_cas9ngg_database \\\n    \
                 --reference chr22.fa.gz --enzyme spcas9ngg\n  \
             java -Xmx4g -jar FlashFry-assembly-1.15.jar discover \\\n    \
                 --database chr22_cas9ngg_database \\\n    \
                 --fasta EMX1_GAGTCCGAGCAGAAGAAGAAGGG.fasta \\\n    \
                 --output EMX1.output\n"
        );
        return;
    }

    let enzyme = Enzyme::spcas9_ngg();
    let finder = SiteFinder::new(enzyme.clone());

    // ---- Step 1: discover guides in EMX1 ----
    let guides_raw = discover_guides(EMX1_PATH, &finder);
    let our_guide_targets: HashSet<String> =
        guides_raw.iter().map(|(s, _)| s.decode_ascii(23)).collect();

    // ---- Step 2: parse FlashFry's discover output ----
    let flashfry_rows = parse_flashfry_output(FLASHFRY_OUTPUT_PATH);
    eprintln!("FlashFry reports {} guides", flashfry_rows.len());
    for row in &flashfry_rows {
        eprintln!(
            "  {} ({}) — otCount={}",
            row.target, row.orientation, row.ot_count
        );
    }

    // ---- Validation A: same guide set ----
    let flashfry_guide_targets: HashSet<String> =
        flashfry_rows.iter().map(|r| r.target.clone()).collect();
    assert_eq!(
        our_guide_targets, flashfry_guide_targets,
        "guide discovery mismatch:\nOurs: {our_guide_targets:?}\nFlashFry: {flashfry_guide_targets:?}",
    );
    eprintln!(
        "✓ guide discovery matches ({} guides)",
        our_guide_targets.len()
    );

    // ---- Validation B: strand orientation matches ----
    for row in &flashfry_rows {
        let expected_strand = match row.orientation.as_str() {
            "FWD" => Strand::Forward,
            "RVS" => Strand::Reverse,
            o => panic!("unexpected FlashFry orientation: {o}"),
        };
        let our_match = guides_raw
            .iter()
            .find(|(s, _)| s.decode_ascii(23) == row.target)
            .expect("we found this guide");
        assert_eq!(
            our_match.1.strand, expected_strand,
            "strand mismatch for {}: FlashFry={}, ours={:?}",
            row.target, row.orientation, our_match.1.strand
        );
    }
    eprintln!("✓ strand orientations match");

    // ---- Step 3: build chr22 BinTable (this is the slow part: ~50Mb FASTA) ----
    let t0 = std::time::Instant::now();
    let (chr22_table, _contig_names) = build_table(CHR22_PATH, enzyme.clone());
    eprintln!(
        "✓ chr22 BinTable built in {:.2}s — {} sites in {} bins",
        t0.elapsed().as_secs_f64(),
        chr22_table.total(),
        chr22_table.occupied_bins(),
    );

    // ---- Step 4: scan chr22 for off-targets of each guide ----
    let guides: Vec<Guide> = guides_raw
        .iter()
        .map(|(s, _)| Guide {
            id: s.decode_ascii(23),
            site: *s,
        })
        .collect();
    let scanner = BinScanner::new(&chr22_table);
    let t1 = std::time::Instant::now();
    let hits = scanner.scan(&guides, MAX_MISMATCHES);
    eprintln!(
        "✓ scan completed in {:.2}s — {} total hits across {} guides",
        t1.elapsed().as_secs_f64(),
        hits.len(),
        guides.len(),
    );

    // ---- Step 5: group hits by guide × unique off-target sequence ----
    // FlashFry collapses identical OT sequences into one entry with a
    // position count, so we compare against the *distinct sequence* count.
    let mut our_distinct_ots: HashMap<usize, HashMap<String, u8>> = HashMap::new();
    for hit in &hits {
        let seq = hit.off_target.decode_ascii(23);
        our_distinct_ots
            .entry(hit.guide_index as usize)
            .or_default()
            .entry(seq)
            .or_insert(hit.mismatches);
    }

    // ---- Validation C: per-guide off-target count + sequence set ----
    for ff_row in &flashfry_rows {
        let our_idx = guides
            .iter()
            .position(|g| g.id == ff_row.target)
            .expect("guide present");
        let ours = our_distinct_ots.get(&our_idx).cloned().unwrap_or_default();
        // FlashFry's `otCount` excludes the on-target (matched 0-mm site).
        // Our scanner returns the on-target as a 0-mm hit; subtract it.
        let our_minus_on_target = ours.iter().filter(|(_, mm)| **mm > 0).count();
        eprintln!(
            "  guide {}: FlashFry ot_count={}, ours(distinct, mm>0)={}",
            ff_row.target, ff_row.ot_count, our_minus_on_target
        );

        assert_eq!(
            our_minus_on_target, ff_row.ot_count,
            "off-target count mismatch for guide {}",
            ff_row.target
        );

        // Also compare the actual sequence sets.
        let our_seqs: HashSet<String> = ours
            .into_iter()
            .filter(|(_, mm)| *mm > 0)
            .map(|(s, _)| s)
            .collect();
        let ff_seqs: HashSet<String> = ff_row
            .off_target_sequences
            .iter()
            .map(|(s, _)| s.clone())
            .collect();
        let only_ours: Vec<&String> = our_seqs.difference(&ff_seqs).collect();
        let only_flashfry: Vec<&String> = ff_seqs.difference(&our_seqs).collect();
        assert!(
            only_ours.is_empty() && only_flashfry.is_empty(),
            "OT sequence set mismatch for {}:\n  only ours: {:?}\n  only FlashFry: {:?}",
            ff_row.target,
            only_ours,
            only_flashfry
        );
    }

    eprintln!("✓ all off-target counts and sequence sets match FlashFry");
}
