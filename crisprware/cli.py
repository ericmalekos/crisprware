#!/usr/bin/env python3
"""CRISPRware command-line interface."""

import argparse
import importlib
import sys

from crisprware.version import __version__


SUBCOMMANDS = {
    "preprocess_annotation": {
        "module_path": "crisprware.preprocess_annotation",
        "help": "Preprocess GTF/GFF annotations with optional RNA-seq filtering.",
    },
    "generate_guides": {
        "module_path": "crisprware.generate_guides",
        "help": "Generate sgRNA sequences matching specified PAM.",
    },
    "index_genome": {
        "module_path": "crisprware.index_genome",
        "help": "Build a crispr-ots off-target index (PAM, protospacer length, and PAM orientation define the enzyme).",
    },
    "score_guides": {
        "module_path": "crisprware.score_guides",
        "help": "Score guides with RS3 cleavage and crispr-ots/Guidescan2 off-target specificity.",
    },
    "rank_guides": {
        "module_path": "crisprware.rank_guides",
        "help": "Rank and select best guides based on scoring criteria.",
    },
}


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="crisprware",
        description="CRISPRware: Tools for CRISPR-based genome editing analysis.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command")

    # Register subcommands with lazy module loading — each module is only
    # imported when its subcommand is actually invoked, so heavy dependencies
    # (pybedtools, rs3, etc.) are not loaded for --version or --help.
    for name, info in SUBCOMMANDS.items():
        subparsers.add_parser(name, help=info["help"], add_help=False)

    # If no arguments, print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # First pass: identify which subcommand was requested
    args, remaining = parser.parse_known_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Import the module and build its full argument parser
    info = SUBCOMMANDS[args.command]
    module = importlib.import_module(info["module_path"])

    sub_parser = argparse.ArgumentParser(
        prog=f"crisprware {args.command}",
        description=info["help"],
    )
    module.add_arguments(sub_parser)
    sub_args = sub_parser.parse_args(remaining)

    module.main(sub_args)


if __name__ == "__main__":
    main()
