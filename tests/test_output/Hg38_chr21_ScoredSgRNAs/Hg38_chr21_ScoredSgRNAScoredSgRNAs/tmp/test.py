#!/usr/bin/env python3

import subprocess

# cmd = [
#     'guidescan',
#     'enumerate',
#     '--max-off-targets',
#     '-1',
#     '--threads',
#     str(threads),
#     '--mismatches',
#     str(mismatches),
#     '--format',
#     'csv',
#     '--rna-bulges',
#     str(rna_bulges),
#     '--dna-bulges',
#     str(dna_bulges),
#     '--threshold',
#     str(threshold),
#     '--mode',
#     mode,
#     "--alt-pam",
#     alt_pam,
#     '--kmers-file',
#     guideCSV,
#     '--output',
#     output,
#     guideIndex
# ]

cmd = [
    'guidescan',
    'enumerate',
    '--threads',
    "8",
    '--format',
    'csv',
    '--threshold',
    '-1',
    '--mode',
    'succint',
    "--alt-pam",
    alt_pam,
    '--kmers-file',
    guideCSV,
    '--output',
    output,
    guideIndex
]

subprocess.run(cmd, check=True)
