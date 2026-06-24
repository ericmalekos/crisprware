# The UCSC Cas12a track pipeline

Four steps take a genome FASTA to a browser-ready track:

1. **Setup** the environment, the GPU engine binary, and the UCSC tools.
2. **Index + guides** build the off-target index (PAM `TTTN`) and enumerate every `TTTV` protospacer.
3. **Score + assemble** run the one `score_guides --ucscgb` command: on-target + off-target scoring,
   internally chunked and fanned across GPUs, merged into one track.
4. **Outputs** inspect the three files the build produces.

```{toctree}
:maxdepth: 1

00_overview
01_setup
02_index_and_guides
03_score_track
04_outputs
```
