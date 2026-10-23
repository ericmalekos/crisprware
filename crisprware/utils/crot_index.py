"""Read the header of a crispr-ots v2 off-target index and check it against an enzyme.

The header is a fixed 128-byte little-endian struct followed by three variable
sections, so reading it costs one 4 KiB read no matter how large the index is.
That is what makes :func:`assert_enzyme_pams` usable as a pre-flight on every
run rather than something a caller might skip for speed.

Why this exists: an AsCas12a result is only correct if the index it was scanned
against contains no TTTT site, because AsCas12a does not cleave a TTTT-prefixed
protospacer. Nothing in the scanner or the track builder could see the index's
PAM set, so three generations of code asserted "TTTV-only index" in comments
while eight genomes were in fact indexed with TTTN. The check below reads the
enumerated PAM list out of the index and refuses the run, which is the one place
that assumption can be turned into a fact.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path

MAGIC = b"CRSPROT2"
HEADER_BYTES = 128

# Field offsets inside DbHeader (crispr-db/src/mmap_db.rs). Named rather than
# unpacked positionally so a future field appended into `_reserved` cannot
# silently shift them.
# 8s magic | 10 x u32 | 10 x u64 (the last two are `_reserved`). The assertion
# below is load-bearing: get the format wrong and the variable-section offsets
# read as zero, which returns an EMPTY PAM list and so silently passes a check
# whose whole job is to reject TTTT.
_H = struct.Struct("<8s10I10Q")
assert _H.size == 128, f"DbHeader format is {_H.size} bytes, expected 128"

#: ``pam_side`` encoding on disk.
PAM_THREE_PRIME = 0
PAM_FIVE_PRIME = 1

#: Enzymes whose off-target set must exclude TTTT, and the PAM each requires.
#: AsCas12a (WT and the 2xNLS variant) cleaves TTTV only; enAsCas12a cleaves
#: TTTT as well, so it is deliberately absent.
TTTV_ONLY_ENZYMES = {"ascas12a"}


@dataclass(frozen=True)
class CrotHeader:
    """The parts of a ``.crot`` header a caller can act on."""

    version: int
    bin_width: int
    pam_side: int
    pam_len: int
    protospacer_len: int
    total_entries: int
    enzyme_name: str
    pams: tuple[str, ...]

    @property
    def pam_side_label(self) -> str:
        return "5'" if self.pam_side == PAM_FIVE_PRIME else "3'"

    def has_tttt(self) -> bool:
        """True if TTTT is among the enumerated PAMs."""
        return any(p.upper() == "TTTT" for p in self.pams)


def resolve_index_path(index: str | Path) -> Path:
    """Accept either the stem crispr-ots is given on the command line or the file itself."""
    p = Path(index)
    return p if p.suffix == ".crot" else p.with_suffix(p.suffix + ".crot")


def read_crot_header(index: str | Path) -> CrotHeader:
    """Read the header of a crispr-ots v2 index.

    Raises:
        FileNotFoundError: no ``.crot`` at ``index``.
        ValueError: the file is not a crispr-ots v2 index.
    """
    path = resolve_index_path(index)
    with open(path, "rb") as fh:
        raw = fh.read(HEADER_BYTES)
        if len(raw) < HEADER_BYTES:
            raise ValueError(f"{path} is shorter than a crispr-ots header ({len(raw)} bytes)")
        (
            magic, version, _flags, bin_width, pam_side, pam_len, protospacer_len,
            pam_count, _contig_count, _bin_count, _pad0,
            total_entries, _mask_lo, _mask_hi, enzyme_name_off, pam_data_off,
            _contig_off, _bin_offsets_off, _entries_off, _res0, _res1,
        ) = _H.unpack_from(raw)
        if magic != MAGIC:
            raise ValueError(f"{path} is not a crispr-ots v2 index (magic {magic!r})")
        fh.seek(enzyme_name_off)
        enzyme_name = fh.read(pam_data_off - enzyme_name_off).decode("ascii", "replace")
        fh.seek(pam_data_off)
        pam_blob = fh.read(pam_count * pam_len).decode("ascii", "replace")
    pams = tuple(pam_blob[i * pam_len:(i + 1) * pam_len] for i in range(pam_count))
    return CrotHeader(
        version=version,
        bin_width=bin_width,
        pam_side=pam_side,
        pam_len=pam_len,
        protospacer_len=protospacer_len,
        total_entries=total_entries,
        enzyme_name=enzyme_name,
        pams=pams,
    )


def assert_enzyme_pams(index: str | Path, enzyme: str) -> CrotHeader:
    """Refuse an index whose enumerated PAMs are wrong for ``enzyme``.

    Returns the header on success so a caller can log or record it.

    Raises:
        ValueError: ``enzyme`` excludes TTTT but the index enumerates it.
    """
    header = read_crot_header(index)
    if enzyme in TTTV_ONLY_ENZYMES and header.has_tttt():
        raise ValueError(
            f"\n\t{resolve_index_path(index)} enumerates TTTT "
            f"(PAMs: {', '.join(header.pams)}; enzyme label '{header.enzyme_name}').\n"
            f"\t{enzyme} does not cleave a TTTT-prefixed protospacer, so every TTTT site in "
            "this index is a false off-target: it inflates the mismatch-bucket counts that "
            "drive the close-off-target grey, the off-target list in the details popup, the "
            "off-target count, max CFD, the saturation cap and the 0-mm screen behind "
            "'unique'.\n"
            "\tRebuild the index with a TTTV PAM:\n"
            f"\t  crisprware index_genome -f <genome.fa> --pam TTTV -l {header.protospacer_len} "
            f"--pam_5_prime --bin_width {header.bin_width}\n"
        )
    return header
