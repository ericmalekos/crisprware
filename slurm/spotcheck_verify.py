#!/usr/bin/env python3
"""Re-derive AsCas12a off-target scoring for a sample of guides, from first principles.

This imports nothing from crisprware. The CFD matrix is parsed from the published
``cas12a_2xnls.csv``, every off-target sequence is pulled out of the genome FASTA with
``samtools faidx``, and the mismatch counts, per-off-target CFD and aggregate
specificity are all recomputed here. Each is then compared against the scanner's own
unfloored output AND against what the finished bigBed track and its crisprDetails popup
publish, so a disagreement anywhere between the index and the browser shows up.

Why it exists: the track asserted a TTTV off-target search for eight genomes while
actually being built on a TTTN index, and nothing in the build or in verify_track.sh
could have caught it, because every check compared the track against itself. This one
compares it against the genome.

Usage (normally via spotcheck_track.sh, which prepares the inputs):
    spotcheck_verify.py <workdir> <genome.fa> <details.tab.gz> <matrix.csv> <cfd_floor> <list_cap>

<workdir> must hold picked.bed (rows sampled from the track), out.csv[.ot.tsv,.guides.tsv]
(an unfloored --cfd-threshold 0 rescan of those guides) and out0.csv (the same guides at
--mismatches 0, which is the only way to see a pre-screened guide's multiplicity).
Exit status is 0 only if every check passes.
"""
import csv, subprocess, sys
from collections import defaultdict

if len(sys.argv) != 7:
    sys.exit(__doc__)
D, FNA, DET, MATRIX, FLOOR, LIST_CAP = sys.argv[1:5] + [float(sys.argv[5]), int(sys.argv[6])]
PROTO = 23
RC = str.maketrans("ACGTN", "TGCAN")
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

# --- matrix -----------------------------------------------------------------------
# penalty[(guide base, off-target protospacer base)][position 0..22], position 1 being
# PAM-proximal. The matrix's d<base> is the base on the strand COMPLEMENTARY to the
# protospacer, so it is complemented on the way in (the Rust does this as `3 - dna` on
# the 2-bit code). Keying on the literal letter instead leaves two thirds of the
# lookups missing, and a missing lookup defaults to a penalty of 1.0, i.e. to "this
# mismatch costs nothing": every CFD comes back as 1.0 and still looks like a number.
pen = defaultdict(lambda: [1.0] * PROTO)
with open(MATRIX) as f:
    r = csv.reader(f)
    next(r)
    for row in r:
        if not row:
            continue
        _, pos, mm, val = row
        g, o = mm.split(":")
        pen[(g[1].replace("U", "T"), COMP[o[1]])][int(pos) - 1] = float(val)
pen = dict(pen)
# A gap would score as 1.0 rather than as an error, so require full coverage.
missing = [(g, o) for g in "ACGT" for o in "ACGT" if g != o and (g, o) not in pen]
if missing:
    sys.exit(f"matrix {MATRIX} is missing {len(missing)} mismatch pairs: {missing}")


def cfd(guide, ot):
    """Product of per-position penalties over the mismatched positions only."""
    s = 1.0
    for i, (g, o) in enumerate(zip(guide, ot)):
        if g != o:
            s *= pen[(g, o)][i]
    return s


def fetch(regions, batch=4000):
    """Sequences for 1-based inclusive 'chrom:start-end' regions, upper-cased."""
    out = []
    for i in range(0, len(regions), batch):
        chunk = regions[i:i + batch]
        p = subprocess.run(["samtools", "faidx", FNA] + chunk,
                           capture_output=True, text=True, check=True)
        seqs, cur = [], []
        for line in p.stdout.splitlines():
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
        if len(seqs) != len(chunk):
            sys.exit(f"faidx returned {len(seqs)} sequences for {len(chunk)} regions")
        out += [s.upper() for s in seqs]
    return out


def popup(offset):
    """The guide's crisprDetails line, by uncompressed byte offset into the bgzip file."""
    if offset == 0:
        return None
    raw = subprocess.run(["bgzip", "-b", str(offset), "-s", "400000", "-c", DET],
                         capture_output=True, check=True).stdout
    return raw.split(b"\n", 1)[0].decode()


# --- inputs -----------------------------------------------------------------------
with open(f"{D}/out.csv.guides.tsv") as f:
    gid2id = {int(r["guide_id"]): r["id"] for r in csv.DictReader(f, delimiter="\t")}
with open(f"{D}/out.csv") as f:
    mode1 = {r["id"]: r for r in csv.DictReader(f)}
with open(f"{D}/out0.csv") as f:
    zero = {r["id"]: r for r in csv.DictReader(f)}

published = {}
with open(f"{D}/picked.bed") as f:
    for p in (l.rstrip("\n").split("\t") for l in f):
        published[f"{p[0]}:{p[1]}:{p[5]}"] = dict(
            score=int(p[4]), guide=p[11], pam=p[12], unique=p[13], spec=p[14],
            offset=int(p[20]))

ots = defaultdict(list)
with open(f"{D}/out.csv.ot.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        ots[gid2id[int(r["guide_id"])]].append(
            (r["chrom"], int(r["start"]), r["strand"], int(r["mm"]), float(r["cfd"])))

# --- per-guide checks -------------------------------------------------------------
bad = 0
n_ot_total = 0
for gid in sorted(published, key=lambda k: (k.rsplit(":", 2)[0], int(k.rsplit(":", 2)[1]))):
    pub, m1 = published[gid], mode1.get(gid)
    guide, rows = pub["guide"], ots[gid]
    print(f"\n=== {gid}  {pub['pam']}|{guide}")
    if m1 is None:
        print("  FAIL: guide absent from the rescan output")
        bad += 1
        continue

    # A --threshold 0 run pre-screens non-unique guides out BEFORE the CFD scan, and a
    # 0-mm hit goes to on_count and returns before the off-target writer, so a dropped
    # guide's multiplicity is only visible as 1/(on_count + 0) in the mm=0 scan.
    if m1["dropped"] == "1":
        z = zero.get(gid)
        n_on = round(1.0 / float(z["specificity_tttv"])) if z else 0
        print(f"  pre-screened out. perfect TTTV matches genome-wide: {n_on} "
              f"(other loci reported: {z['mismatch_counts'] if z else '-'})")
        if n_on < 2:
            print("  FAIL: dropped but fewer than 2 perfect matches")
            bad += 1
        else:
            print("  PASS: >1 perfect TTTV match, so the 0-mm screen is right to drop it")
        if pub["unique"] != "False" or pub["spec"] != "" or pub["offset"] != 0:
            print(f"  FAIL: expected unique=False / blank spec / _offset 0, got "
                  f"{pub['unique']} / {pub['spec']!r} / {pub['offset']}")
            bad += 1
        else:
            print("  PASS track row: unique=False, blank specificity, _offset 0")
        continue

    seqs = fetch([f"{c}:{s + 1}-{s + 27}" for c, s, *_ in rows]) if rows else []
    pams, protos = [], []
    for (c, s, strand, mm, cf), seq in zip(rows, seqs):
        if strand == "-":
            seq = seq.translate(RC)[::-1]
        pams.append(seq[:4])
        protos.append(seq[4:4 + PROTO])
    n_ot_total += len(rows)

    # 1. TTTV. The point of the whole exercise: no TTTT may appear anywhere.
    counts = {p: pams.count(p) for p in ("TTTA", "TTTC", "TTTG", "TTTT")}
    offenders = [p for p in pams if not (p[:3] == "TTT" and p[3] in "ACG")]
    print(f"  off-targets reported: {len(rows)}   PAMs: "
          + " ".join(f"{k}={v}" for k, v in counts.items() if v))
    if offenders:
        print(f"  FAIL: {counts['TTTT']} TTTT and {len(offenders) - counts['TTTT']} "
              f"other non-TTTV PAMs")
        bad += 1
    else:
        print("  PASS TTTV: every reported off-target sits behind TTTA/TTTC/TTTG")

    # 2. mismatch count and CFD, re-derived from the genome.
    mm_bad = cfd_bad = 0
    worst = 0.0
    for (c, s, st, mm, cf), proto in zip(rows, protos):
        if sum(1 for a, b in zip(guide, proto) if a != b) != mm:
            mm_bad += 1
        my = cfd(guide, proto)
        if abs(my - cf) > 2e-6 * max(1.0, cf):
            cfd_bad += 1
        worst = max(worst, abs(my - cf))
    if rows:
        print(f"  {'PASS' if not mm_bad else 'FAIL'} mismatches: "
              f"{len(rows) - mm_bad}/{len(rows)} agree")
        print(f"  {'PASS' if not cfd_bad else 'FAIL'} CFD: "
              f"{len(rows) - cfd_bad}/{len(rows)} agree (max |delta| {worst:.2e})")
        bad += bool(mm_bad) + bool(cfd_bad)

    # 3. the specificity arithmetic: 1/(on_count + sum over mm>0), guidescan convention.
    on = sum(1 for r in rows if r[3] == 0) or 1
    ssum = sum(cfd(guide, p) for (c, s, st, mm, cf), p in zip(rows, protos) if mm > 0)
    my_spec = 1.0 / (on + ssum)
    print(f"  on-target loci: {on}   sum(CFD over mm>0): {ssum:.6f}")
    print(f"  re-derived 1/({on} + sum) = {my_spec:.6f}   "
          f"scanner says {m1['specificity_tttv']}")
    if abs(my_spec - float(m1["specificity_tttv"])) > 1e-5:
        print("  FAIL: specificity does not reproduce")
        bad += 1
    else:
        print("  PASS specificity")
    if m1["specificity_tttn"] != m1["specificity_tttv"]:
        print(f"  FAIL: tttn {m1['specificity_tttn']} != tttv, so the index holds TTTT")
        bad += 1

    # 4. what the track published.
    print(f"  track: unique={pub['unique']}, spec={pub['spec'] or '-'}, score={pub['score']}")
    if pub["unique"] != "True":
        print("  FAIL: unique flag disagrees with the rescan")
        bad += 1
    elif abs(round(my_spec, 3) - float(pub["spec"])) > 1.1e-3:
        print(f"  FAIL: published spec {pub['spec']} != round({my_spec:.6f}, 3)")
        bad += 1
    else:
        print("  PASS published specificity matches to 3 dp")

    # 5. the popup list and its histogram, both floored at the build's CFD threshold.
    line = popup(pub["offset"])
    if line is None:
        print(f"  popup: none (_offset 0), off_target_count={m1['off_target_count']}")
        if int(m1["off_target_count"] or 0) != 0:
            print("  FAIL: _offset 0 but the guide has off-targets")
            bad += 1
        continue
    hist, lst = line.split("\t")
    nbins = len(hist.split(","))
    keep = [(c, s, st, cf) for (c, s, st, mm, cf) in rows if mm > 0 and cf >= FLOOR]
    mine = [f"{c};{s}{st};{int(round(cf * 1000))}" for c, s, st, cf in keep]
    theirs = lst.split("|") if lst else []
    myhist = [0] * nbins
    for (c, s, st, mm, cf) in rows:
        if mm > 0 and cf >= FLOOR and mm < nbins:
            myhist[mm] += 1
    print(f"  popup: {len(theirs)} listed, histogram {hist}")
    if len(mine) > LIST_CAP:
        # Over list_cap the build emits an empty list on purpose and keeps the true
        # histogram, so the browser shows its "too many off-targets" message.
        if theirs:
            print(f"  FAIL: {len(mine)} off-targets is over list_cap {LIST_CAP} "
                  f"but {len(theirs)} were listed")
            bad += 1
        else:
            print(f"  PASS over list_cap {LIST_CAP}: list blank, histogram kept")
    elif set(mine) != set(theirs):
        only_m, only_t = set(mine) - set(theirs), set(theirs) - set(mine)
        print(f"  FAIL popup list: {len(only_m)} missing, {len(only_t)} extra")
        for x in list(only_m)[:3] + list(only_t)[:3]:
            print("    ", x)
        bad += 1
    else:
        print("  PASS popup list: same off-targets, same rounded CFDs")
    if ",".join(map(str, myhist)) != hist:
        print(f"  FAIL histogram: re-derived {','.join(map(str, myhist))}")
        bad += 1
    else:
        print("  PASS histogram matches the floored list")

print(f"\n==== {len(published)} guides, {n_ot_total:,} off-targets re-derived from the "
      f"genome: {'ALL CHECKS PASSED' if bad == 0 else str(bad) + ' FAILURES'} ====")
sys.exit(1 if bad else 0)
