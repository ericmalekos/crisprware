#!/usr/bin/env python3

"""
    This script takes the output from generate_guides and adds
    with  on-target Ruleset3 and off-target Guidescan2 scores.

    score_guides --grna_bed tests/test_data/gRNA_test.bed \
        --guidescan2_indices tests/test_data/Hg38_chr21_gscan_index/chr21.fa.index\
        --tracr Chen2013 --threads 2 --output_prefix tests/test_output/Hg38_chr21

"""

from typing import List, Optional, Union

import pandas as pd
import subprocess
import argparse
from rs3.seq import predict_seq

# from scipy.stats import norm
from pathlib import Path
from os import remove
from crisprware.utils.utility_functions import create_output
from crisprware.utils.dna_sequence_functions import map_ambiguous_sequence


def restricted_float(x: Union[str, float]) -> float:
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return x


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add score_guides arguments to the given parser."""
    parser.add_argument("-b", "--grna_bed", type=str, help="grnas.bed ouput of generate_guides.", required=True)

    parser.add_argument(
        "-i",
        "--guidescan2_indices",
        type=str,
        help="One or more, space-separate Guidescan2 indices. \
            A specificity score will be calculated against each index separately.",
        required=False,
        nargs="*",
    )

    parser.add_argument(
        "--tracr",
        type=str,
        default=None,
        required=False,
        choices=["Hsu2013", "Chen2013", "both"],
        help="TracrRNA version for cleavage scoring. \
            Either 'Hsu2013' or 'Chen2013' or 'both', see https://github.com/gpp-rnd/rs3 \
            for details.",
    )

    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads [default: 8]")

    parser.add_argument(
        "--threshold",
        type=int,
        default=2,
        help="Threshold for Guidescan2 off-target hits. If off-targets are found this distance away \
         the gRNA will be discarded, i.e. set to 2 to discard any guides with a 0, 1 or 2 mismatches \
         from another PAM adjacent sequence. --threshold=-1 to retain all guides [default: 2]",
    )

    parser.add_argument(
        "--mismatches", type=int, default=3, help="Number of mismatches for Guidescan2 off-target scoring [default: 3]"
    )

    parser.add_argument(
        "--rna_bulges", type=int, default=0, help="RNA bulges for Guidescan2 off-target scoring [default: 0]"
    )

    parser.add_argument(
        "--dna_bulges", type=int, default=0, help="DNA bulges for Guidescan2 off-target scoring [default: 0]"
    )

    parser.add_argument(
        "--mode",
        type=str,
        default="succinct",
        choices=["succinct", "complete"],
        help="Whether Guidescan2 temporary output should be succinct or complete mode [default: 0]",
    )

    parser.add_argument(
        "--alt_pams",
        type=str,
        help="One or more, space-separate alternative pams for off-target \
            consideration. e.g. NAG",
        required=False,
        nargs="*",
    )

    parser.add_argument(
        "-d",
        "--drop_duplicates",
        help="Drop exact duplicate gRNAs before scoring to save time.\
             Set flag to retain duplicates. [default: True]",
        default=True,
        action="store_false",
    )

    parser.add_argument(
        "--skip_rs3", help="Set flag to skip RS3 scoring [default: False]", default=False, action="store_true"
    )

    parser.add_argument(
        "--skip_gs2", help="Set flag to skip Guidescan2 scoring [default: False]", default=False, action="store_true"
    )

    parser.add_argument(
        "--min_rs3",
        type=float,
        default=float("-inf"),
        help="Minimum cleavage RS3 score. RS3 cleavage scores are formatted \
            as z-scores, so this is interpreted as a standard deviation cutoff. \
            Functionality also available in rank_guides.py. Applying at this \
            stage can increase speed by filtering before off-target scoring. \
            [default: None]",
    )

    parser.add_argument(
        "--cas12a_scorer",
        choices=["none", "enpam_gb", "deepcpf1", "enseq_deepcpf1", "seq_deepcpf1variants", "both"],
        default="none",
        help="Cas12a on-target scorer. enpam_gb for en(As)Cas12a; deepcpf1 \
            for wildtype AsCas12a/LbCas12a (Kim 2018); enseq_deepcpf1 for \
            wildtype AsCas12a (Chen 2025, modern); seq_deepcpf1variants for \
            variant-specific scoring (requires --cas12a_variant). \
            'both' runs enpam_gb + deepcpf1 in parallel. Mutually exclusive \
            with --tracr (RS3 is SpCas9-only); implies --skip_rs3. [default: none]",
    )

    parser.add_argument(
        "--cas12a_variant",
        type=str,
        default=None,
        help="Cas12a variant name (e.g. AsCas12a_Ultra, enAsCas12a-HF1, \
            LbCas12a, HyperFi-AsCas12a). Required when --cas12a_scorer is \
            seq_deepcpf1variants. See crisprware.scorers.seq_deepcpf1variants \
            for the full 23-variant list.",
    )

    parser.add_argument(
        "--min_deepcpf1",
        type=float,
        default=float("-inf"),
        help="Minimum DeepCpf1 score (raw regression, ~[0, 100]). Applied \
            after scoring; analogous to --min_rs3. [default: None]",
    )

    parser.add_argument(
        "--min_enpam_gb",
        type=float,
        default=float("-inf"),
        help="Minimum enPAM+GB score (probability-like, [0, 1]). Applied \
            after scoring; analogous to --min_rs3. [default: None]",
    )

    parser.add_argument(
        "--min_enseq_deepcpf1",
        type=float,
        default=float("-inf"),
        help="Minimum enseq-DeepCpf1 / seq-DeepCpf1variants score \
            (probability, [0, 1]). Applied after scoring. [default: None]",
    )

    parser.add_argument(
        "--cas9_scorer",
        choices=["none", "deepspcas9", "deephf_wt_u6", "deephf_esp", "deephf_hf"],
        default="none",
        help="Additional SpCas9 on-target scorer to run alongside RS3. \
            deepspcas9: Kim 2019 inception-CNN (Sci Adv); 30-nt context \
            (4 + 20 protospacer + 3 PAM + 3 downstream); unbounded \
            regression. deephf_*: Wang 2019 BiLSTM (Nat Commun) for three \
            Cas9 variants -- wildtype SpCas9 (wt_u6), eSpCas9 (esp), \
            SpCas9-HF1 (hf); 23-nt protospacer+PAM input; output in [0,1]. \
            Runs in parallel with --tracr (no mutex). [default: none]",
    )

    parser.add_argument(
        "--min_deepspcas9",
        type=float,
        default=float("-inf"),
        help="Minimum DeepSpCas9 score (unbounded regression, ~[0, 100]). \
            Applied after scoring; analogous to --min_rs3. [default: None]",
    )

    parser.add_argument(
        "--min_deephf",
        type=float,
        default=float("-inf"),
        help="Minimum DeepHF score (probability, [0, 1]). Applied to whichever \
            deephf_* variant is selected. [default: None]",
    )

    parser.add_argument(
        "--chunk_size",
        type=int,
        default=100000,
        help="Number of gRNAs to hold in memory for cleavage scoring and off-target filtering. \
        Reduce if memory constrained. Increasing may improve runtime [default: 100000]",
    )

    parser.add_argument(
        "-k",
        "--keep_tmp",
        help="Set flag to keep temporary Guidescan2 output [default: False]",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "-o", "--output_directory", help="Path to output. [default: current directory]", type=str, default=""
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scores gRNAs from generate_guides.")
    add_arguments(parser)
    return parser.parse_args()


def gscan_scoring(
    guideCSV: str,
    output: str,
    guideIndex: str,
    threads: int = 2,
    mismatches: int = 3,
    rna_bulges: int = 0,
    dna_bulges: int = 0,
    threshold: int = 2,
    mode: str = "succinct",
    alt_pam: str = "",
    keep_gscan: bool = False,
) -> pd.DataFrame:
    """
    Executes GuideScan's enumerate command with specified parameters and processes the output.

    Parameters:
    - guideCSV (str): The path to the input CSV file containing guide RNA sequences.
    - output (str): The path where the output CSV file will be saved.
    - guideIndex (str): The index used for off-target scoring.
    - threads (int, optional): The number of threads to use for the GuideScan command. Defaults to 2.
    - mismatches (int, optional): The maximum number of mismatches allowed. Defaults to 3.
    - rna_bulges (int, optional): The number of RNA bulges allowed. Defaults to 0.
    - dna_bulges (int, optional): The number of DNA bulges allowed. Defaults to 0.
    - threshold (int, optional): The threshold parameter for removing low mismatches. Defaults to 2.
    - mode (str, optional): The mode of operation for the GuideScan command. Defaults to "succinct".
    - alt_pam (str, optional): An alternative PAM sequence to be considered. Defaults to an empty string.

    Returns:
    - pandas.DataFrame: A DataFrame containing the processed output from the GuideScan command.
    """

    cmd = [
        "crispr-ots",
        "enumerate",
        "--max-off-targets",
        "-1",
        "--threads",
        str(threads),
        "--mismatches",
        str(mismatches),
        "--format",
        "csv",
        "--rna-bulges",
        str(rna_bulges),
        "--dna-bulges",
        str(dna_bulges),
        "--threshold",
        str(threshold),
        "--mode",
        mode,
        "--alt-pam",
        "None",
        "--kmers-file",
        guideCSV,
        "--output",
        output,
        guideIndex,
    ]

    subprocess.run(cmd, check=True)

    # read the csv file
    gscanDF = pd.read_csv(output)

    # gscanDF = process_gscan(gscanDF, guideIndex)
    if gscanDF.empty:
        print('\n\tGuidescan2 output is empty, consider lowering "--threshold"\n')

    gscanDF = gscanDF.drop_duplicates(subset="id", keep="first")

    gscanDF["specificity"] = gscanDF["specificity"].round(4)

    db_name = Path(guideIndex).parts[-1]

    gscanDF.rename(columns={"specificity": "specificity_" + db_name}, inplace=True)

    gscanDF = gscanDF.drop(["sequence", "match_chrm", "match_position", "match_strand", "match_distance"], axis=1)

    if not keep_gscan:
        remove(output)
        remove(guideCSV)

    return gscanDF


def compute_rs3_scores(gRNAlist: List[str], tracr: str, threads: int, chunk_size: int) -> List[float]:
    """
    Compute RS3 scores for a list of gRNAs with the specified tracrRNA sequence.

    Parameters:
    - gRNAlist (list): List of gRNA 30-nt contexts.
    - tracr (str): The tracrRNA sequence identifier to be used in cleavage scoring.
    - threads (int): The number of parallel jobs to run for cleavage scoring.
    - chunk_size (int): The size of chunks to split the gRNA list into for processing.

    Returns:
    - list: List of RS3 scores for the input gRNA sequences.
    """
    gRNAScores = []

    # Iterate over big_list in chunks of size chunk_size
    for i in range(0, len(gRNAlist), chunk_size):
        sublist = gRNAlist[i : i + chunk_size]
        processed_sublist = predict_seq(sublist, sequence_tracr=tracr, n_jobs=threads)
        gRNAScores.extend(processed_sublist)

    return gRNAScores


def cleavage_scoring(
    gRNADF: pd.DataFrame, tracr: str, threads: int = 2, chunk_size: int = 200000, minStdDev: float = float("-inf")
) -> pd.DataFrame:
    """
    Computes RS3 cleavage scores for gRNAs using the specified tracrRNA sequence and filters based on minimum standard deviation.

    Parameters:
    - gRNADF (pandas.DataFrame): DataFrame containing gRNA sequences with a column named 'context' for gRNA 30-nt contexts.
    - tracr (str): The tracrRNA sequence identifier to be used in cleavage scoring.
    - threads (int, optional): The number of parallel jobs to run for cleavage scoring. Defaults to 2.
    - chunk_size (int, optional): The size of chunks to split the gRNA list into for processing, to manage memory usage. Defaults to 2000000.
    - minStdDev (float or None, optional): The minimum standard deviation threshold for filtering gRNAs based on their RS3 z-score.
      gRNAs with a z-score below this threshold will be excluded. Defaults to None, which disables filtering.

    Returns:
    - pandas.DataFrame: The input DataFrame augmented with 'RS3_score' and 'rs3_cdf' columns, and optionally filtered
      based on the 'RS3_score' threshold.
    """

    print("\n\tBeginning RS3 cleavage scoring\n\tIf memory constrained reduce '--chunk_size'\n")

    gRNAlist = gRNADF["context"].tolist()

    if tracr == "both":
        gRNAScores_Hsu2013 = compute_rs3_scores(gRNAlist, "Hsu2013", threads, chunk_size)
        gRNAScores_Chen2013 = compute_rs3_scores(gRNAlist, "Chen2013", threads, chunk_size)

        gRNADF["RS3_score_Hsu2013"] = gRNAScores_Hsu2013
        gRNADF["RS3_score_Hsu2013"] = gRNADF["RS3_score_Hsu2013"].round(4)
        # gRNADF['rs3_cdf_Hsu2013'] = norm.cdf(gRNADF['RS3_score_Hsu2013'])
        # gRNADF['rs3_cdf_Hsu2013'] = gRNADF['rs3_cdf_Hsu2013'].round(4)

        gRNADF["RS3_score_Chen2013"] = gRNAScores_Chen2013
        gRNADF["RS3_score_Chen2013"] = gRNADF["RS3_score_Chen2013"].round(4)
        # gRNADF['rs3_cdf_Chen2013'] = norm.cdf(gRNADF['RS3_score_Chen2013'])
        # gRNADF['rs3_cdf_Chen2013'] = gRNADF['rs3_cdf_Chen2013'].round(4)
    elif tracr in ["Hsu2013", "Chen2013"]:
        gRNAScores = compute_rs3_scores(gRNAlist, tracr, threads, chunk_size)

        gRNADF["RS3_score_" + tracr] = gRNAScores
        gRNADF["RS3_score_" + tracr] = gRNADF["RS3_score_" + tracr].round(4)
        # gRNADF['rs3_cdf'] = norm.cdf(gRNADF['RS3_score'])
        # gRNADF['rs3_cdf'] = gRNADF['rs3_cdf'].round(4)
        gRNADF = gRNADF[gRNADF["RS3_score_" + tracr] > minStdDev]

    gRNADF = gRNADF.copy()

    return gRNADF


def check_files_exist(index: str) -> None:
    """Check the existence of the three required files for a given index."""

    files = [f"{index}.reverse", f"{index}.forward", f"{index}.gs"]
    for file in files:
        if not Path(file).exists():
            raise FileNotFoundError(
                f"\n\n\tRequired file {file} not found for index {index}. \
                                    \n\tMake sure \n\t\t{index}.reverse \n\t\t{index}.forward \n\t\t{index}.gs \n\texist\n"
            )


def get_alt_pams(pams: List[str]) -> str:
    """
    Generates a string of alternative PAM sequences by expanding each ambiguous PAM sequence in the input list.

    Parameters:
    - pams (list of str): A list of PAM sequences. These sequences can include ambiguous nucleotide symbols
      that represent multiple possible nucleotides at a particular position.

    Returns:
    - str: A string containing all unique alternative PAM sequences derived from the input list, space-separated.
    """

    pamlist = []
    for apam in pams:
        pamlist += map_ambiguous_sequence(apam)
    pamstr = " ".join(set(pamlist))
    return pamstr


def main(args: Optional[argparse.Namespace] = None) -> None:

    if args is None:
        args = parse_arguments()

    # print("ALT PAMS" + args.alt_pams)
    if not args.skip_gs2:
        for index in args.guidescan2_indices:
            check_files_exist(index)

    if args.cas12a_scorer != "none":
        if args.tracr:
            raise ValueError(
                "\n\t--tracr (SpCas9 RuleSet3) and --cas12a_scorer are mutually exclusive."
                " For Cas12a guides, drop --tracr.\n"
            )
        args.skip_rs3 = True

    if args.cas12a_scorer == "seq_deepcpf1variants" and not args.cas12a_variant:
        raise ValueError("\n\t--cas12a_variant is required when --cas12a_scorer is seq_deepcpf1variants.\n")

    if not args.skip_rs3 and not args.tracr:
        raise ValueError(
            "\n\tEither --tracr (SpCas9 RS3) or --cas12a_scorer must be set, or pass --skip_rs3 to skip on-target scoring entirely.\n"
        )

    pams = ""
    if args.alt_pams:
        pams = get_alt_pams(args.alt_pams)
        print(f"\n\tConsidering {pams} for off-target scoring \n")
    else:
        args.alt_pams = None

    # gRNA_output_path = "./" + args.output_prefix + "ScoredSgRNAs/" + args.output_prefix.split("/")[-1] + "ScoredSgRNAs.tsv"
    # tmp_path = create_output_directory(base_dir="./" + args.output_prefix + "ScoredSgRNAs/",output_prefix="tmp/")

    gRNA_output_path, tmp_path = create_output(
        args.grna_bed, outdir=args.output_directory, extension="scoredgRNA", stripped="_gRNA", tmp=True
    )
    gRNA_output_path += ".bed"
    # print(gRNA_output_path)
    # print(tmp_path)

    specificity_cols = []

    gRNADF = pd.read_csv(args.grna_bed, delimiter="\t", header=0)
    final_columns = gRNADF.columns.tolist() + ["sequence"]
    final_columns.remove("id,sequence,pam,chromosome,position,sense")

    if args.drop_duplicates:
        print(f"\n\tBefore dropping duplicates:\t{len(gRNADF)}")
        gRNADF["gRNA"] = gRNADF.iloc[:, 3].str.split(",").str[1]

        # Step 3 & 4: Find unique and duplicate k-mers, then filter the DataFrame to keep only the unique ones
        unique_mask = ~gRNADF["gRNA"].duplicated(keep=False)
        gRNADF = gRNADF[unique_mask]
        print(f"\tAfter dropping duplicates:\t{len(gRNADF)}\n")

    if not args.skip_rs3:
        if args.tracr == "both":
            final_columns += ["RS3_score_Hsu2013", "RS3_score_Chen2013"]
        else:
            final_columns += ["RS3_score_" + args.tracr]
        gRNADF = cleavage_scoring(
            gRNADF=gRNADF, tracr=args.tracr, chunk_size=args.chunk_size, threads=args.threads, minStdDev=args.min_rs3
        )
        print(f"\n\tAfter dropping RS3 cleavage scores below {args.min_rs3}:\t{len(gRNADF)}\n")

    if args.cas12a_scorer in ("enpam_gb", "both"):
        from crisprware.scorers import enpam_gb as _enpam_gb

        print("\n\tBeginning enPAM+GB Cas12a on-target scoring\n")
        gRNADF["enpam_gb_score"] = _enpam_gb.compute_enpam_gb_scores(
            gRNADF["context"].tolist(), threads=args.threads, chunk_size=args.chunk_size
        )
        gRNADF["enpam_gb_score"] = gRNADF["enpam_gb_score"].round(8)
        if args.min_enpam_gb > float("-inf"):
            before = len(gRNADF)
            gRNADF = gRNADF[gRNADF["enpam_gb_score"] > args.min_enpam_gb]
            print(f"\tAfter dropping enPAM+GB scores below {args.min_enpam_gb}: {before} -> {len(gRNADF)}\n")
        final_columns += ["enpam_gb_score"]

    if args.cas12a_scorer in ("deepcpf1", "both"):
        from crisprware.scorers import deepcpf1 as _deepcpf1

        print("\n\tBeginning DeepCpf1 Cas12a on-target scoring\n")
        gRNADF["deepcpf1_score"] = _deepcpf1.compute_deepcpf1_scores(
            gRNADF["context"].tolist(), threads=args.threads, chunk_size=args.chunk_size
        )
        gRNADF["deepcpf1_score"] = gRNADF["deepcpf1_score"].round(8)
        if args.min_deepcpf1 > float("-inf"):
            before = len(gRNADF)
            gRNADF = gRNADF[gRNADF["deepcpf1_score"] > args.min_deepcpf1]
            print(f"\tAfter dropping DeepCpf1 scores below {args.min_deepcpf1}: {before} -> {len(gRNADF)}\n")
        final_columns += ["deepcpf1_score"]

    if args.cas12a_scorer == "enseq_deepcpf1":
        from crisprware.scorers import enseq_deepcpf1 as _enseq

        print("\n\tBeginning enseq-DeepCpf1 Cas12a on-target scoring\n")
        gRNADF["enseq_deepcpf1_score"] = _enseq.compute_enseq_deepcpf1_scores(
            gRNADF["context"].tolist(), threads=args.threads, chunk_size=args.chunk_size
        )
        gRNADF["enseq_deepcpf1_score"] = gRNADF["enseq_deepcpf1_score"].round(8)
        if args.min_enseq_deepcpf1 > float("-inf"):
            before = len(gRNADF)
            gRNADF = gRNADF[gRNADF["enseq_deepcpf1_score"] > args.min_enseq_deepcpf1]
            print(
                f"\tAfter dropping enseq-DeepCpf1 scores below {args.min_enseq_deepcpf1}: {before} -> {len(gRNADF)}\n"
            )
        final_columns += ["enseq_deepcpf1_score"]

    if args.cas12a_scorer == "seq_deepcpf1variants":
        from crisprware.scorers import seq_deepcpf1variants as _variants

        canonical = _variants._normalize_variant(args.cas12a_variant)
        col_name = f"seq_deepcpf1variants_{canonical.replace('-', '_')}_score"
        print(f"\n\tBeginning seq-DeepCpf1variants scoring (variant={canonical})\n")
        gRNADF[col_name] = _variants.compute_variant_scores(
            gRNADF["context"].tolist(),
            variant=canonical,
            threads=args.threads,
            chunk_size=args.chunk_size,
        )
        gRNADF[col_name] = gRNADF[col_name].round(8)
        if args.min_enseq_deepcpf1 > float("-inf"):
            before = len(gRNADF)
            gRNADF = gRNADF[gRNADF[col_name] > args.min_enseq_deepcpf1]
            print(f"\tAfter dropping scores below {args.min_enseq_deepcpf1}: {before} -> {len(gRNADF)}\n")
        final_columns += [col_name]

    if args.cas9_scorer == "deepspcas9":
        from crisprware.scorers import deepspcas9 as _deepspcas9

        print("\n\tBeginning DeepSpCas9 SpCas9 on-target scoring\n")
        gRNADF["deepspcas9_score"] = _deepspcas9.compute_deepspcas9_scores(
            gRNADF["context"].tolist(), threads=args.threads, chunk_size=args.chunk_size
        )
        gRNADF["deepspcas9_score"] = gRNADF["deepspcas9_score"].round(8)
        if args.min_deepspcas9 > float("-inf"):
            before = len(gRNADF)
            gRNADF = gRNADF[gRNADF["deepspcas9_score"] > args.min_deepspcas9]
            print(f"\tAfter dropping DeepSpCas9 scores below {args.min_deepspcas9}: {before} -> {len(gRNADF)}\n")
        final_columns += ["deepspcas9_score"]

    if args.cas9_scorer.startswith("deephf_"):
        from crisprware.scorers import deephf as _deephf

        variant = args.cas9_scorer.removeprefix("deephf_")  # wt_u6, esp, or hf
        col_name = f"deephf_{variant}_score"
        print(f"\n\tBeginning DeepHF SpCas9 on-target scoring (variant={variant})\n")
        # DeepHF takes 23-nt protospacer+PAM. The pam column in the id field
        # is generic ("NGG"); the actual 3-bp PAM bases live in the context
        # column immediately after the 20-bp protospacer. Slice 23 nt out by
        # finding the protospacer (from the sequence column) inside context.
        composite = gRNADF["id,sequence,pam,chromosome,position,sense"].str.split(",")
        protospacers = composite.str[1]

        def _extract_23(row):
            pro = row["protospacer"]
            ctx = row["context"]
            i = ctx.find(pro)
            if i < 0 or i + 23 > len(ctx):
                return None
            return ctx[i:i + 23]

        tmp_df = pd.DataFrame({"protospacer": protospacers.tolist(), "context": gRNADF["context"].tolist()})
        seq23 = tmp_df.apply(_extract_23, axis=1)
        gRNADF[col_name] = _deephf.compute_deephf_scores(
            seq23.tolist(), variant=variant, threads=args.threads, chunk_size=args.chunk_size
        )
        gRNADF[col_name] = gRNADF[col_name].round(8)
        if args.min_deephf > float("-inf"):
            before = len(gRNADF)
            gRNADF = gRNADF[gRNADF[col_name] > args.min_deephf]
            print(f"\tAfter dropping DeepHF scores below {args.min_deephf}: {before} -> {len(gRNADF)}\n")
        final_columns += [col_name]

    gRNADF.loc[:, "id"] = gRNADF["id,sequence,pam,chromosome,position,sense"].str.split(",").str[0]
    guidescan_dfs = []

    if not args.skip_gs2:
        for gscanIndex in args.guidescan2_indices:
            print(
                "\n\tBeginning Guidescan2 specificity scoring against "
                + gscanIndex
                + "\n\tIf memory constrained reduce '--chunk_size'\n"
            )
            guidescan_chunk_dfs = []
            suffix_index = gscanIndex.split("/")[-1]
            # chunk up the input and save it in a form compatible with guidescan processing
            for i, (_, chunk) in enumerate(gRNADF.groupby(gRNADF.index // args.chunk_size)):
                input = tmp_path + f"{suffix_index}Input.{i + 1}.csv"
                print("input:" + input)
                output = input.replace("Input", "Output")
                print(f"\n\tSaved Guidescan input file to {input}\n")
                chunk[["id,sequence,pam,chromosome,position,sense"]].to_csv(input, sep="\t", index=False)
                guidescan_chunk_dfs.append(
                    gscan_scoring(
                        guideCSV=input,
                        output=output,
                        guideIndex=gscanIndex,
                        threads=args.threads,
                        mismatches=args.mismatches,
                        rna_bulges=args.rna_bulges,
                        dna_bulges=args.dna_bulges,
                        mode=args.mode,
                        threshold=args.threshold,
                        alt_pam=pams,
                        keep_gscan=args.keep_tmp,
                    )
                )

            guidescan_dfs.append(pd.concat(guidescan_chunk_dfs, ignore_index=True))

        for df in guidescan_dfs:
            gRNADF = gRNADF.merge(df, on="id", how="outer")
            # gRNADF = gRNADF.merge(df, on='id', how='inner')

        # If any rows have NaN for all specificity columns, drop those rows
        specificity_cols = [col for col in gRNADF.columns if "specificity" in col]
        rows_to_drop = gRNADF[specificity_cols].isna().all(axis=1)
        gRNADF = gRNADF.loc[~rows_to_drop]

    gRNADF[["sequence"]] = gRNADF["id,sequence,pam,chromosome,position,sense"].str.split(",", expand=True).iloc[:, 1:2]

    gRNADF = gRNADF.drop(["id,sequence,pam,chromosome,position,sense", "id"], axis=1)

    gRNADF = gRNADF[final_columns + specificity_cols]

    gRNADF.to_csv(gRNA_output_path, na_rep="-1", sep="\t", index=False)
    print(f"\n\tSaved output file to {gRNA_output_path}\n")


if __name__ == "__main__":
    main()
