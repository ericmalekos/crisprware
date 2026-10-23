#!/bin/bash
# Assert the structural invariants of a crispr-ots aggregated (Mode-1) output.
#
# One streaming awk pass, so it runs on a multi-GB agg. Works on both layouts, detected
# from the header rather than from the filename:
#
#   SpCas9   id,specificity,max_cfd,off_target_count,mismatch_counts,saturated,dropped
#   Cas12a   id,specificity_tttn,specificity_tttv,max_cfd,off_target_count,
#            mismatch_counts,saturated,dropped
#
# Checks, in order of what they would have caught:
#
#   1. mismatch_counts is POPULATED whenever off_target_count > 0. This is the check that
#      SpCas9 under --scanner gpu failed silently: every other column was right and the
#      histogram was 0;0;0;0, so nothing errored and the column read as "no off-targets at
#      any mismatch level" for 1.46M of 1.50M rows.
#   2. sum(mismatch_counts[1..]) == off_target_count. Both come off the same accumulator and
#      a cap truncates them together, so this holds for saturated rows too.
#   3. mismatch_counts[0] == 0. The mm-0 bucket is emitted from on_count, not from the
#      histogram; a nonzero value here means a 0-mm hit leaked into the off-target path.
#   4. Cas12a only: specificity_tttn == specificity_tttv on every scored row. These can
#      differ ONLY if the index enumerates a TTTT PAM, so the equality is a cheap
#      index-level assertion that the off-target space really was TTTV. Eight published
#      genomes would have failed this.
#
# Not checked here: whether the CFD values are right. That needs the genome, and
# spotcheck_track.sh does it by re-deriving every off-target from the FASTA.
#
# Usage:  assert_agg_invariants.sh <agg.csv> [more.csv ...]
# Exit 0 only if every file passes. --quiet prints one line per file.
set -uo pipefail

QUIET=0
[ "${1:-}" = "--quiet" ] && { QUIET=1; shift; }
[ $# -gt 0 ] || { echo "usage: assert_agg_invariants.sh [--quiet] <agg.csv> [...]" >&2; exit 2; }

rc=0
for f in "$@"; do
  [ -s "$f" ] || { echo "FATAL: missing or empty $f" >&2; rc=1; continue; }
  awk -F, -v FILE="$f" -v QUIET="$QUIET" '
    NR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      # A Cas12a agg is the one carrying both specificity flavours.
      cas12a = ("specificity_tttn" in col) && ("specificity_tttv" in col)
      need = "off_target_count mismatch_counts saturated dropped"
      n = split(need, want, " ")
      for (i = 1; i <= n; i++)
        if (!(want[i] in col)) { printf "FATAL: %s has no %s column\n", FILE, want[i]; bad_hdr = 1 }
      if (bad_hdr) exit 3
      OTC = col["off_target_count"]; MMC = col["mismatch_counts"]
      SAT = col["saturated"];        DRP = col["dropped"]
      if (cas12a) { TTN = col["specificity_tttn"]; TTV = col["specificity_tttv"] }
      next
    }
    {
      rows++
      if ($DRP == 1) { dropped++; next }
      scored++
      if ($SAT == 1) saturated++
      otc = $OTC + 0

      # 1. populated whenever there is anything to report
      if (otc > 0 && ($MMC == "" || $MMC ~ /^0(;0)*$/)) { unpop++; if (!ex1) ex1 = $1 }

      # 2 and 3. buckets against the total, and the unused mm-0 bucket
      if ($MMC != "") {
        nb = split($MMC, b, ";")
        s = 0
        for (i = 2; i <= nb; i++) s += b[i] + 0
        if (s != otc) { sumbad++; if (!ex2) ex2 = sprintf("%s otc=%d counts=%s", $1, otc, $MMC) }
        if (b[1] + 0 != 0) { zerobad++; if (!ex3) ex3 = sprintf("%s counts=%s", $1, $MMC) }
      }

      # 4. the TTTV index assertion
      if (cas12a && $TTN != "" && $TTV != "" && $TTN != $TTV) {
        tttbad++; if (!ex4) ex4 = sprintf("%s tttn=%s tttv=%s", $1, $TTN, $TTV)
      }
    }
    END {
      fail = unpop + sumbad + zerobad + tttbad
      if (QUIET == 1) {
        printf "%s  rows=%d scored=%d dropped=%d saturated=%d  %s\n", FILE, rows, scored,
               dropped, saturated, (fail ? "FAIL" : "ok")
      } else {
        printf "\n=== %s%s\n", FILE, (cas12a ? "  (Cas12a, two specificity flavours)" : "  (SpCas9)")
        printf "  rows %d | scored %d | dropped %d | saturated %d\n", rows, scored, dropped, saturated
        printf "  %s unpopulated mismatch_counts with off_target_count > 0: %d\n",
               (unpop ? "FAIL" : "PASS"), unpop
        if (unpop) printf "       e.g. %s\n", ex1
        printf "  %s sum(buckets[1..]) != off_target_count: %d\n", (sumbad ? "FAIL" : "PASS"), sumbad
        if (sumbad) printf "       e.g. %s\n", ex2
        printf "  %s nonzero mm-0 bucket: %d\n", (zerobad ? "FAIL" : "PASS"), zerobad
        if (zerobad) printf "       e.g. %s\n", ex3
        if (cas12a) {
          printf "  %s specificity_tttn != specificity_tttv (index holds TTTT): %d\n",
                 (tttbad ? "FAIL" : "PASS"), tttbad
          if (tttbad) printf "       e.g. %s\n", ex4
        }
      }
      exit (fail ? 1 : 0)
    }
  ' "$f" || rc=1
done

[ "$QUIET" = 1 ] || { echo; [ $rc -eq 0 ] && echo "=== ALL FILES PASS ===" || echo "=== FAILURES ==="; }
exit $rc
