# Project notes — caveats from FlashFry and GuideScan2

Things learned during exploration that affect design decisions. Keep this
file current as new gotchas surface during implementation.

## FlashFry caveats

**Hard 24-base cap on encoded sites.** Single `Long`, 48 bits sequence + 16
bits count. `bitcoding/BitEncoding.scala:47, 202`. Any guide+PAM > 24 bp
overflows. This is why FlashFry can't do 23-nt Cas12a (4+23 = 27).

**Enzymes are hard-coded `case object`s.** `standards/StandardScanParameters.scala`
defines six fixed parameter packs (`Cas9ParameterPack`, `Cpf1ParameterPack`,
…). PAM regex is literal-letter alternation (e.g. `[ACGTacgt]{20}[AG]G`), not
IUPAC. Adding a PAM requires a code change + recompile.

**Comparison masks are hard-coded per enzyme.**
`Cas9ParameterPack.comparisonBitEncoding = 0x3FFFFFFFFFC0L` etc. — these
literally pin 20-nt guide + 3-bp PAM. Any length change needs a new mask
constant.

**Bit-encoding direction.** First character of the string ends up in the
*highest* used bits (because encoding does `<<2` before OR'ing each base).
For 3' PAM (Cas9): low 6 bits = PAM. For 5' PAM (Cpf1): high 8 bits within
the 48-bit window = PAM. The two layouts are *not* interchangeable; bin
slicing in `crispr/BinWriter.scala:59-66` keys off this.

**CFD scoring is per-off-target string decode.**
`scoring/Doench2016CFDScore.scala:60-77, 132-150`: decodes Long → String for
every off-target, walks `.zip`, does `HashMap[String, Double]` lookup with
keys like `"rC:dG,9"`. Major hot-path hazard. Replace with a flat 4×4×L
lookup indexed by 2-bit guide/OT base and position.

**BGZF + ByteBuffer per block.** `reference/traverser/LinearTraverser.scala:122-130`,
`utils/Utils.scala:154-186`. Each block read allocates a fresh `byte[]`,
decompresses via HTSJDK's BGZF, then copies into a `Long[]` via
`ByteBuffer.getLong()` in a loop. For NVMe-class storage this dominates
runtime.

**Bin lookup hashes a String.** `BinaryHeader.blockOffsets: HashMap[String, BlockOffset]`.
Bin key is a 6- or 7-character String hashed every access. Use the bin's
encoded `u32` as a direct array index.

**Single-threaded by design.** No `.par` collections, no futures in the
scan path. The only "parallelism" is the overflow callback culling.

**No bulge support.** FlashFry only does Hamming distance. Cas9 tolerates
small indels in vivo; if that matters, FlashFry is the wrong baseline.

## GuideScan2 caveats

**`--start` flag is run-global, not per-guide.** A single boolean
(`opts.start`, default false = 3' PAM) for the entire run. Cannot mix Cas9
and Cas12a guides in one query. `guidescan.cxx:48,53`,
`process.hpp:63,66,71,81`, `printer.hpp:194,212,325`.

**CFD is hard-wired to 3' PAM.** `printer.hpp:142` does
`pam = match_sequence.substr(20, 3)` — pulls PAM from positions 20-22 of
the match. If you run with `--start` (5' PAM), specificity is
undefined/broken. The `mm_scores` and `pam_scores` tables in `doench.hpp`
are SpCas9-only by construction.

**CFD also assumes 20-nt protospacer.** `printer.hpp:99`:
`if ((sgRNA.length() != 20) || (PAM.length() != 3)) return 1.0` — returns
unity (no penalty) on any other length. So GuideScan2 effectively can't
score Cas12a even though it can *find* its off-targets.

**Per-recursive-call heap allocation.** `index.hpp:147,166,226,244` —
`match + c` / `sequence + c` does a `std::string` concat per recursion node.
With a 20-nt query and k=3 mismatches this fans to thousands of allocs per
guide. A fixed-size scratch buffer eliminates it entirely.

**`std::function` callback dispatch.** `index.hpp:121,132,191,262` — every
recursive call goes through type-erased `std::function`, blocking inlining
of the callback. Templated functor would be free.

**`fmt::format` in CFD inner loop.** `printer.hpp:107` builds a key string
per mismatch, then looks it up in `unordered_map<string, float>`. Same
hazard FlashFry has, different language.

**`std::set<match>` accumulator.** Log-time inserts + RB-tree allocs per
match. `std::vector` + sort/unique at the end is faster.

**Two FM-indexes (fwd + revcomp), not a bidirectional one.** Doubles index
size on disk. Bidirectional BWT exists in SDSL but isn't used.

**Index is PAM-agnostic but k-mer extraction is a separate Python script.**
`scripts/generate_kmers.py` does the FASTA-scan + PAM-expansion (TTTN →
{TTTA,TTTC,TTTG,TTTT}) and emits a CSV of candidate sites. The C++ binary
only consumes that CSV; it doesn't scan FASTAs itself. This is a flexibility
win (any PAM works) but it means PAM logic lives in two languages.

## Shared caveats (apply to both)

**CFD is a SpCas9-specific model.** The Doench 2016 paper trained CFD on
SpCas9 with NGG. No equivalent published table exists for Cas12a, xCas9,
SpRY, etc. For non-SpCas9 enzymes our tool should either skip CFD or note
that the score is extrapolated and meaningless. Don't promise CFD on
Cas12a.

**Output positions are 1-indexed in both.** Plan on 1-indexed in our CLI
output for FlashFry parity and SAM convention; keep internal coords
0-indexed.

**Strand encoding differs.** FlashFry uses `F`/`R` in TSV output; GuideScan2
uses `+`/`-` in CSV and SAM FLAG 0/16. We should default to `+`/`-` for SAM
compatibility.

## Science caveats

**CFD ≠ MIT/Hsu.** Different mismatch models, different score scales, not
interchangeable. Users may expect one or the other; document which we emit.

**Specificity is a derived score, not raw CFD.** Both tools compute
`1 / (1 + Σ CFD_off_target × count_at_position)`. The numerator/denominator
choice differs slightly between papers; transcribe FlashFry's formula
exactly to maintain output equivalence.

**FlashFry vs GuideScan2 specificity disagree by exactly the on-target
multiplicity.** FlashFry/ours: `1 / (1 + Σ cfd × count)` with the on-
target as a singular `+1` in the denominator. GuideScan2 emits one
CFD-1.0 row per *on-target position*, so its denominator is
`Σ_all_rows(cfd) ≈ n_on_target + Σ_off_target(cfd)`. The two reconcile
exactly: `spec_gs2 = 1 / (1/spec_ours + (n_on_target - 1))`. Verified
end-to-end by `cli_guidescan2_equivalence` on the 1000-gRNA chr22
fixture; max |Δ predicted - actual| = 6e-7 (GuideScan2's printed
precision floor). We follow FlashFry's convention by default.

**GuideScan2's default PAM matching is exact, not NGG.** Without
`--alt-pam`, `guidescan enumerate` only matches off-targets whose
PAM equals the guide's `pam` column in the input CSV byte-for-byte.
For NGG-equivalent search, pass `--alt-pam AGG --alt-pam CGG
--alt-pam GGG --alt-pam TGG`. Forgetting this causes catastrophic
spec disagreement (max |Δ| ~0.9 on the 1000-gRNA fixture).

**Off-target counts at k=4+ explode.** Most real workflows cap at k=3.
Be willing to refuse or warn at high k. FlashFry has a max-off-target
overflow mechanism (`crispr/CRISPRSiteOT.scala`); replicate it.

## Reproducible 1000-gRNA fixture generation

```python
# Generate /tmp/ffry-quickstart/random_1000.fasta and
# /tmp/ffry-quickstart/random_1000.gs2.csv from chr22.fa.gz.
import gzip, re, random
random.seed(42)

with gzip.open('/tmp/ffry-quickstart/chr22.fa.gz', 'rt') as f:
    seq = ''.join(l.strip() for l in f if not l.startswith('>')).upper()

starts = [m.start() for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', seq)]
sampled = sorted(random.sample(starts, 1000))

with open('/tmp/ffry-quickstart/random_1000.fasta', 'w') as out:
    for i, start in enumerate(sampled):
        out.write(f'>chr22_site_{i}\n{seq[start:start+23]}\n')

with open('/tmp/ffry-quickstart/random_1000.gs2.csv', 'w') as out:
    out.write('id,sequence,pam,chromosome,position,sense\n')
    for i, start in enumerate(sampled):
        site = seq[start:start+23]
        out.write(f'chr22_site_{i},{site[:20]},{site[20:]},chr22,{start+1},+\n')
```

## Open questions to revisit during implementation

- How does FlashFry handle ambiguous bases (Ns) in the genome? Skipped at
  index time or encoded as a fifth symbol? Verify in
  `reference/ReferenceEncoder.scala`.
- Does FlashFry's `OrderedBinTraversalFactory` actually save time vs. the
  linear traversal in practice? Profile before reimplementing the smart
  traversal.
- Threshold for switching between linear-block and indexed-block storage in
  FlashFry's DB is `maxTargetsPerLinearBin = 500` (`DatabaseWriter.scala:85`).
  Validate that 500 is still the right cutoff on modern hardware, or pick
  fresh from a microbenchmark.
