# Off-target records

A guide's `_offset` is the uncompressed byte position of its line in the track's `*_crisprDetails.tab`
(`ascas12a_crisprDetails.tab` / `encas12a_crisprDetails.tab`, named for the
track), resolved through the `.gzi` index. Each line is two columns:

```text
_mismatchCounts   _crisprOfftargets
0,2,0,0,2         chr2;182255206-;466|chr15;101948968-;344|chr19;83578+;344|chr10;98762096+;172
```

- **`_mismatchCounts`** off-target counts by mismatch bucket `[0,1,2,3,4]` over the 23-nt guide; the
  on-target site is not counted. Raw absolute totals, unaffected by chunking.
- **`_crisprOfftargets`** `|`-separated. Each entry is `chrom;pos±;cfd`, a single CFD score x1000 from
  that track's matrix; every track scans exactly one matrix. An off-target is listed
  only if its CFD clears `--ucscgb_cfd_threshold`; the list is capped at `--ucscgb_list_cap`. (So counts
  can exceed the number listed.) On the AsCas12a track a `TTTT`-PAM off-target has its AsCas12a CFD zeroed
  (AsCas12a does not cleave it), so it can appear with a `0` score.

## Coordinate geometry

`pos` is the 0-based **leftmost** base of the 27-nt window (PAM + spacer); `±` is the strand the guide
binds. The browser fetches `[pos, pos+27)`, reverse-complements it on `-`, then reads the first 4 nt as
the PAM and the next 23 as the match, the same layout as the on-target row.

## `_offset = 0` and uniqueness

`_offset = 0` means **no off-target line is stored** and the browser shows *"No off-targets found."* Two
distinct cases produce it:

| Case | `unique` | Why no list |
|------|----------|-------------|
| Non-unique guide | `False` | an identical protospacer exists elsewhere; off-target analysis is skipped |
| Unique, perfectly specific | `True` | specificity `1.0`, nothing cleared the CFD floor |

So **uniqueness is read from the `unique` flag, not from `_offset`.** Example non-unique row:
`TTTATTTTTTTTTTGTAGAGATG`, `unique = False`, `_offset = 0`.
