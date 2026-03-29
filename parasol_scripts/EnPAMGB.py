#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EnPAM_GB scoring utility.

Reads a TSV of candidate sgRNAs, scores the sequences in a specified column
(default: "context") using the EnPAM_GB model from `sgrna_modeler`, and writes
an output TSV with a new score column.

Typical usage
-------------
python EnPAMGB.py \
  -i sgRNAs/sgRNA.6.tsv \
  -o sgRNAs/sgRNA.6.scored.tsv \
  --seq-col context \
  --out-col enpam_gb_score

Input format
------------
- A tab-delimited file with at least one column containing sequences to score.
- If the first header cell begins with '#', e.g. "#chr", it will be normalized
  (leading '#' stripped) for convenience.

Output
------
- The same table as input, plus a new column with EnPAM_GB scores.
"""

import argparse
from typing import Iterable, List

import pandas as pd
from sgrna_modeler import models as sg
from sgrna_modeler import enzymes as en


def getEnPAMGB(sequences: Iterable[str]) -> List[float]:
    """
    Score a list of sequences with the EnPAM_GB model (Cas12a).

    Parameters
    ----------
    sequences : Iterable[str]
        Iterable of nucleotide sequences to score. Each entry should be a
        string representing the context sequence expected by the model.

    Returns
    -------
    List[float]
        A list of scores in the same order as `sequences`.

    Notes
    -----
    - Internally constructs a `SklearnSgrnaModel`, loads EnPAM_GB weights via
      `sg.get_enpam_gb()`, and predicts with `model.predict_seqs`.
    - This function initializes a model for each call; for large inputs or
      repeated calls, prefer loading the model once (see `score_file`).
    """
    model = sg.SklearnSgrnaModel()
    model_weights = sg.get_enpam_gb()
    model.load_model(model_weights, en.cas12a, "enPAM_GB")
    return list(model.predict_seqs(sequences))


def score_file(
    input_path: str,
    output_path: str,
    seq_col: str = "context",
    out_col: str = "enpam_gb_score",
    batch_size: int = 4096,
) -> None:
    """
    Read a TSV, score sequences from `seq_col` with EnPAM_GB, and write results.

    Parameters
    ----------
    input_path : str
        Path to the input TSV file. Must contain a header row.
    output_path : str
        Path to the output TSV file to write. Will overwrite if it exists.
    seq_col : str, optional
        Column name that contains sequences to score (default: "context").
    out_col : str, optional
        Name of the new column to store EnPAM_GB scores (default: "enpam_gb_score").
    batch_size : int, optional
        Number of sequences to process per batch (default: 4096). Adjust for
        memory/performance trade-offs.

    Returns
    -------
    None
        Writes a TSV to `output_path`.

    Raises
    ------
    ValueError
        If `seq_col` is not present in the input table.

    Notes
    -----
    - If the first header cell begins with '#', it is normalized by stripping
      the leading '#'. This is common in genomic TSVs (e.g., "#chr" -> "chr").
    - Sequences are cast to strings with NaNs replaced by empty strings.
    - The model is instantiated and loaded once per file for efficiency.
    """
    # Read TSV as strings to preserve content verbatim
    df = pd.read_csv(input_path, sep="\t", dtype=str)

    # Normalize first header if it starts with "#"
    if len(df.columns) > 0 and isinstance(df.columns[0], str) and df.columns[0].startswith("#"):
        df.rename(columns={df.columns[0]: df.columns[0].lstrip("#")}, inplace=True)

    if seq_col not in df.columns:
        raise ValueError(f"Column '{seq_col}' not found. Available: {list(df.columns)}")

    # Prepare sequences
    seqs = df[seq_col].fillna("").astype(str).tolist()

    # Load model once
    model = sg.SklearnSgrnaModel()
    model_weights = sg.get_enpam_gb()
    model.load_model(model_weights, en.cas12a, "enPAM_GB")

    # Predict in batches
    scores: List[float] = []
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i : i + batch_size]
        scores.extend(model.predict_seqs(batch))

    # Attach and write
    df[out_col] = scores
    
	# Round to 8 decimal places
    df[out_col] = [round(float(s), 8) for s in scores]
    
    df.to_csv(output_path, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes:
        - input (str):  Input TSV path.
        - output (str): Output TSV path.
        - seq_col (str): Column containing sequences to score.
        - out_col (str): Name of the output score column.
        - batch_size (int): Batch size for scoring.
    """
    parser = argparse.ArgumentParser(
        description="Score sgRNA context sequences with EnPAM_GB and append results."
    )
    parser.add_argument("-i", "--input", required=True, help="Input TSV path")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path")
    parser.add_argument(
        "--seq-col", default="context", help="Column with sequences to score (default: context)"
    )
    parser.add_argument(
        "--out-col", default="enpam_gb_score", help="Name of output score column (default: enpam_gb_score)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=4096, help="Batch size for scoring (default: 4096)"
    )
    return parser.parse_args()


def main() -> None:
    """
    Entry point for CLI execution.

    Reads CLI args, scores the input table, and writes the output TSV.
    Prints the output path upon success.
    """
    args = parse_args()
    score_file(args.input, args.output, args.seq_col, args.out_col, args.batch_size)
    print(f"Wrote scored TSV to {args.output}")


if __name__ == "__main__":
    main()
