## Jutils for Tool-Agnostic Output

**Goal:** Visualize differential splicing output uniformly across rMATS, leafcutter, MntJULiP, and MAJIQ.

**Approach:** Convert tool output to Jutils' standard TSV, then plot. Run from the Jutils clone (`python3 jutils.py ...`). Run on rMATS output; the leafcutter/MntJULiP/MAJIQ converters follow `jutils.py convert-results --help` and were not run.

```bash
scripts/jutils_pipeline.sh /path/to/Jutils rmats_output/ meta.tsv bam_list.tsv annotation.gtf chr1:1000-2000 jutils_out/
```

`meta.tsv` is `sample<TAB>condition` and `bam_list.tsv` is `sample<TAB>bam<TAB>condition`. The script runs, in order: `convert-results` (writes `rmats_JC_results.tsv` and `rmats_JCEC_results.tsv`), `heatmap` (`--q-value 0.05`, env `Q`; **needs >= 2 events passing the cutoffs**, otherwise it prints "Skipping"; writes `clustermap*.pdf`), `sashimi` for the coordinate, and `venn-diagram`. `--tsv-file-list` is a FILE with one `path<TAB>label` line per TSV, not a comma-separated list; the script writes it. Outputs go to `hm/`, `sh/`, `vn/` under the output directory.

(Yang 2021 *Bioinformatics*) Useful when comparing multiple tools' outputs across publications or doing meta-analysis. The sashimi labels are per-sample junction counts and matched pysam.
