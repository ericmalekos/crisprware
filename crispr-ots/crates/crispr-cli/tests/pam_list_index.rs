//! End-to-end round-trip for an explicit multi-PAM index.
//!
//! The `.crot` format stores the enzyme's fully-expanded PAM set, so an enzyme
//! built from a PAM list (whose accepted PAMs are not a single IUPAC pattern)
//! must survive a build/persist/open cycle unchanged. This is the path exercised
//! by `crispr-ots index --pam-list ...`, used for the enAsCas12a tier 1+2+3 set
//! of 47 four-mers.

use std::collections::BTreeSet;

use crispr_cli::discover::enzyme_from_pam_list;
use crispr_cli::{run_build, BuildConfig};
use crispr_db::MmapDb;

fn pam_strings(db: &MmapDb) -> BTreeSet<String> {
    db.enzyme()
        .canonical_pams
        .iter()
        .map(|p| p.iter().map(|b| b.to_ascii() as char).collect())
        .collect()
}

#[test]
fn multi_pam_set_round_trips_through_the_crot() {
    // A Cas12a-geometry set that is deliberately not one IUPAC pattern:
    // TTTA (canonical), TTCC and TATC (tier-1 expanded PAMs).
    let pams = ["TTTA".to_string(), "TTCC".to_string(), "TATC".to_string()];
    let enzyme = enzyme_from_pam_list(&pams, 20, true).expect("valid PAM list");

    // A small reference containing one site per PAM on the forward strand:
    // each is PAM(4) + 20 nt protospacer. Distinct spacers so nothing collides.
    let fasta = "\
>chr_test
TTTAACACACACACACACACACACGGGGGG\
TTCCGTGTGTGTGTGTGTGTGTGTCCCCCC\
TATCTGTGTGTGTGTGTGTGTGTGAAAAAA
";
    let dir = std::env::temp_dir().join(format!("crispr-pamlist-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, fasta).expect("write fasta");
    let crot = dir.join("idx.crot");

    run_build(&BuildConfig {
        reference: fa,
        enzyme: enzyme.clone(),
        output: crot.clone(),
        bin_width: None,
    })
    .expect("build multi-PAM index");

    let db = MmapDb::open(&crot).expect("open .crot");
    let read = db.enzyme();

    // Metadata round-trips exactly.
    assert_eq!(read.pam_len, 4);
    assert_eq!(read.protospacer_len, 20);
    assert_eq!(read.pam_side, enzyme.pam_side);
    assert_eq!(read.compare_mask, enzyme.compare_mask);
    assert_eq!(
        pam_strings(&db),
        BTreeSet::from([
            "TTTA".to_string(),
            "TTCC".to_string(),
            "TATC".to_string()
        ]),
        "the full PAM set must survive the build/open cycle"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn single_and_list_indexes_of_the_same_pams_agree() {
    // Building {TTTA,TTTC,TTTG,TTTT} via --pam-list must yield the same PAM set
    // as the IUPAC TTTN, proving the list path is a strict generalisation.
    let via_list = enzyme_from_pam_list(
        &[
            "TTTA".to_string(),
            "TTTC".to_string(),
            "TTTG".to_string(),
            "TTTT".to_string(),
        ],
        20,
        true,
    )
    .expect("list parses");
    let via_iupac =
        crispr_cli::discover::enzyme_from_pam("TTTN", 20, true).expect("IUPAC parses");

    let set = |e: &crispr_encoding::Enzyme| -> BTreeSet<String> {
        e.canonical_pams
            .iter()
            .map(|p| p.iter().map(|b| b.to_ascii() as char).collect())
            .collect()
    };
    assert_eq!(set(&via_list), set(&via_iupac));
    assert_eq!(via_list.compare_mask, via_iupac.compare_mask);
}
