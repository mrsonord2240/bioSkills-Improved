---
name: bio-sashimi-plots
description: Creates sashimi-style plots showing RNA-seq read coverage and splice junction counts using ggsashimi (general-purpose, condition-grouped overlays), rmats2sashimiplot (rMATS-output-aware), MAJIQ-VOILA (LSV posteriors, interactive viewer; licence-gated), leafviz (leafcutter clusters Shiny), Jutils (tool-agnostic heatmaps and sashimi for rMATS/leafcutter/MntJULiP/MAJIQ output), or pyGenomeTracks (multi-track publication figures). Tool choice depends on the upstream differential-splicing tool's output format and the publication vs interactive use case. Use when visualizing specific splicing events, validating differential splicing calls, or producing publication-quality figures.
tool_type: python
primary_tool: ggsashimi
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with (2026-09-20): ggsashimi 1.1.5, rmats2sashimiplot 4.0.0, leafcutter/leafviz 0.2.9, pyGenomeTracks 3.9, Jutils 1.5, pysam 0.24, pandas 2.2+, R 4.2.3 with **ggplot2 3.4.4**. MAJIQ/VOILA is licence-gated and was not installed: its commands below follow MAJIQ's public docs and were not run.

**ggsashimi needs ggplot2 < 3.5.** With ggplot2 3.5.2 and 4.0.3 the gene-model track and the per-group panels are shifted against each other and x tick labels are clipped (arc counts stay correct); with 3.4.4 everything aligns. **Always look at the rendered figure** for exon edges lining up with arcs before reporting it.

Install (ggsashimi and Jutils are on neither conda nor PyPI; PyPI `jutils` is an unrelated package):

```bash
conda install -c conda-forge -c bioconda rmats2sashimiplot pygenometracks pysam bedtools regtools samtools
conda install -c conda-forge r-base=4.2 r-ggplot2=3.4.4 r-data.table r-gridextra r-gtable seaborn scikit-learn
git clone https://github.com/guigolab/ggsashimi     # ggsashimi.py: one script, needs pysam and R on PATH
export PATH=$PWD/ggsashimi:$PATH                    # the recipes call `ggsashimi.py` by name; its `#!/usr/bin/env python` must find a python with pysam
git clone https://github.com/splicebox/Jutils       # jutils.py
```

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Sashimi Plot Visualization

Visualize RNA-seq coverage tracks with splice junction arcs labeled by read count. Sashimi plots originated with MISO (Katz 2010 *Nat Methods*); modern tools differ in input handling, group aggregation logic, and customization. Tool choice is not interchangeable — some tools work only with specific upstream output formats.

## Tool Selection Matrix

| Tool | Best for | Input | Strengths | Fails when |
|------|----------|-------|-----------|------------|
| ggsashimi | Publication-quality grouped overlays from any BAM | BAMs + region | `--overlay` aggregates samples within a group; clean PDFs | No native rMATS/MAJIQ integration; need to extract coords manually; ggplot2 >= 3.5 misaligns panels |
| rmats2sashimiplot | One-line plot from rMATS output | rMATS event file + BAMs | No manual coord extraction | rMATS-specific; doesn't handle leafcutter or MAJIQ |
| MAJIQ-VOILA | Interactive LSV browsing with posterior PSI distributions | MAJIQ build + psi/deltapsi | Splice-graph topology; LSV-aware; posterior violins | Static figures; licence-gated download, not run in testing |
| leafviz | Cluster-level interactive browsing | leafcutter differential output | Filter table + sashimi-like plots | leafcutter-specific |
| Jutils | Unified output across rMATS, leafcutter, MntJULiP, MAJIQ | Tool-specific differential output | Heatmaps, Venn, sashimi tool-agnostically | Output less polished than ggsashimi |
| pyGenomeTracks | Multi-track publication figures (RNA-seq + ChIP/ATAC) | bedGraph or BigWig + BED + GTF | Combine RNA with chromatin tracks | Not splicing-specific; configure tracks manually |
| IGV (interactive) | Quick ad-hoc inspection | BAM + region | Scrollable, instant | Not for publication figures |
| MISO sashimi | Historical | MISO output | Original sashimi format | MISO unmaintained; no longer recommended |

## Decision Tree by Goal

| Goal | Recommended tool |
|------|-------------------|
| Validate a specific rMATS hit | rmats2sashimiplot (one-line; `references/rmats2sashimiplot.md`) or ggsashimi (custom) |
| Validate a leafcutter cluster | leafviz (interactive; `references/leafviz.md`) or ggsashimi with cluster coordinates |
| Validate a MAJIQ LSV (complex topology) | MAJIQ-VOILA (only tool that shows full LSV graph; `references/majiq-voila.md`) |
| Publication-quality two-condition comparison | ggsashimi `-O 3 -A mean_j` for grouped overlay |
| Multi-track figure (RNA-seq + H3K4me3 + ATAC) | pyGenomeTracks (`references/pygenometracks.md`) |
| Quick ad-hoc browsing during development | IGV sashimi |
| Tool-agnostic batch heatmap of significant events | Jutils (`references/jutils.md`) |
| Interactive cohort-level filtering of leafcutter results | leafviz Shiny (`references/leafviz.md`) |

## ggsashimi for Publication Overlays

**Goal:** Generate publication-quality sashimi plot for a region with samples grouped by condition and per-sample tracks aggregated.

**Approach:** Define samples + groups in a TSV (no header), a palette file, then call ggsashimi with coordinates, GTF, and visual flags. ggsashimi exits 0 on several failures (see Silent Failures), so check the inputs before and the figure after.

```python
import subprocess
import pandas as pd
from pathlib import Path

# ggsashimi input: col1 = sample id, col2 = BAM path, col3 = group (used by -O overlay and -C colour)
groups = pd.DataFrame({
    'sample_id': ['ctrl1', 'ctrl2', 'ctrl3', 'trt1', 'trt2', 'trt3'],
    'bam': ['ctrl1.bam', 'ctrl2.bam', 'ctrl3.bam', 'trt1.bam', 'trt2.bam', 'trt3.bam'],
    'group': ['Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment']
})
groups.to_csv('sashimi_groups.tsv', sep='\t', index=False, header=False)
Path('palette.txt').write_text('#1f77b4\n#ff7f0e\n')  # one colour per group, in order of first appearance

missing = [b for b in groups['bam'] if not Path(b).is_file()]
assert not missing, f'ggsashimi would drop these BAMs silently: {missing}'

subprocess.run([
    'ggsashimi.py',
    '-b', 'sashimi_groups.tsv',
    '-c', 'chr17:43094000-43125000',   # contig spelled as in the BAM header
    '-o', 'BRCA1_sashimi',
    '--alpha', '0.25',
    '--height', '3',
    '--width', '10',
    '--shrink',
    '--fix-y-scale',
    '--ann-height', '4',
    '-g', 'gencode_v45.gtf',
    '--base-size', '14',
    '-O', '3', '-C', '3', '-P', 'palette.txt',
    '-A', 'mean_j',
    '-F', 'pdf'
], check=True)
assert Path('BRCA1_sashimi.pdf').is_file() and Path('BRCA1_sashimi.pdf').stat().st_size > 0, 'no figure written (R error above?)'
```

`examples/plot_sashimi.py` wraps this as `plot_sashimi()` and adds contig mapping (chrX vs X), an empty-region check and the `--shrink` guard below. Checked: ggsashimi's junction labels equal an independent pysam count (planted 3v3 set, real ENCODE 12-BAM locus, real chrX BAMs).

Key ggsashimi flags (Garrido-Martin 2018 *PLoS Comput Biol*):
- `-O 3`: column 3 of the TSV is the overlay level; samples of a group are drawn in one track. Required for `-A`
- `-C 3 -P palette.txt`: colour by column 3 using the palette file (R colour names or hex, one per line). Without `-C` everything is grey; with `-C` and no `-P` the colours are R defaults (red/green), not blue/orange
- `-A mean_j`: the arc label is the rounded (half to even: 6.5 shows 6) plain mean of the raw junction counts of the group's samples that have the junction; every sample's coverage stays overlaid. `mean` and `median` also aggregate the coverage; `median_j` is the median analogue. There is no depth normalization
- `-M`: minimum reads for a junction to be drawn, **default 1**, inclusive, applied **per sample before `-A`**: samples below `-M` drop out of the mean, biasing labels upward (planted 10.33 shows 11 at `-M 10`; 10 at `-M 1`). Keep `-M 1` for exact labels; raise it (5-10, 20+ for crowded plots) only to declutter. If `-M` is above every junction's count, ggsashimi still exits 0 and writes a coverage-only figure with no arcs (`examples/plot_sashimi.py` warns)
- `--shrink`: rescale long introns for compact display (keeps exons visible in genes with multi-kb introns). Crashes with `RuntimeError: generator raised StopIteration` when no junction passes `-M` (with `-s`, when one strand has none): lower `-M` or drop `--shrink`
- `--fix-y-scale`: identical y-axis across groups (essential for visual comparison)
- `--alpha 0.25`: transparency of per-sample coverage in overlay mode
- `--height`/`--width`/`--ann-height`/`--base-size`: sizes in inches / font size (e.g. `--width 12 --height 4 --base-size 16`)
- `-F pdf|svg|png|jpeg|tiff` (`-R` for raster PPI); pick explicitly: PDF for publication, PNG for slides, SVG for editing
- `-g`: GTF with exons; ggsashimi has no feature filter, so pre-filter the GTF (`awk '$0 ~ /protein_coding/'`) to restrict it
- Colour convention: control = blue (`#1f77b4`), treatment = orange (`#ff7f0e`); document it, use ColorBrewer for >2 groups

## Batch Plotting from rMATS Hits

**Goal:** Auto-generate sashimi plots for all significant rMATS differential events.

**Approach:** Continues the ggsashimi block above (it reuses `groups`, `sashimi_groups.tsv` and `palette.txt`). Parse SE.MATS.JC.txt, expand coordinates to flanking exons + 500nt context, map the contig name onto the BAM header, iterate ggsashimi and check every figure. rMATS writes `chrX`; an Ensembl-style BAM calls it `X` and ggsashimi dies with `ValueError: invalid contig`. Events near a contig start give a start < 1.

```python
import re
import subprocess
import pandas as pd
import pysam
from pathlib import Path

contigs = set(pysam.AlignmentFile(groups['bam'][0]).references)  # groups = the TSV above

def bam_contig(name):
    for cand in (name, name.removeprefix('chr'), 'chr' + name.removeprefix('chr')):
        if cand in contigs:
            return cand
    raise ValueError(f'contig {name} not in the BAM header')

diff = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')
sig = diff[(diff['FDR'] < 0.05) & (diff['IncLevelDifference'].abs() > 0.10)]

Path('sashimi_plots').mkdir(exist_ok=True)
failed = []
for _, ev in sig.head(25).iterrows():
    region = f'{bam_contig(ev["chr"])}:{max(1, ev["upstreamES"] - 500)}-{ev["downstreamEE"] + 500}'
    safe_name = re.sub(r'[^A-Za-z0-9._-]', '_', f'{ev["geneSymbol"]}_{ev["chr"]}_{ev["upstreamES"]}_{ev["ID"]}')
    out = Path(f'sashimi_plots/{safe_name}.pdf')
    subprocess.run([
        'ggsashimi.py', '-b', 'sashimi_groups.tsv', '-c', region,
        '-o', str(out.with_suffix('')), '-M', '1', '--shrink', '--fix-y-scale',
        '-O', '3', '-C', '3', '-P', 'palette.txt', '-A', 'mean_j', '-g', 'annotation.gtf', '-F', 'pdf'
    ])
    if not (out.is_file() and out.stat().st_size > 0):
        failed.append(region)
assert not failed, f'no figure for {failed}'
```

MXE files also carry `upstreamES`/`downstreamEE`, and that span already covers both alternative exons. `examples/plot_sashimi.py` `batch_plot_rmats_events()` is the same recipe with the `--shrink` guard and a `RuntimeError` listing every failed event.

## Reference Files

The ggsashimi recipes, interpretation guide, failure modes and Common Errors stay in this file. Read the tool's file when the request needs that tool:

| File | Read when |
|------|-----------|
| `references/rmats2sashimiplot.md` | Plotting rMATS events directly, `--group-info`, `--event-type`, or its exit-0 failures |
| `references/majiq-voila.md` | Browsing MAJIQ LSVs with `voila view` (V2 and V3 inputs; licence-gated, not run) |
| `references/leafviz.md` | The leafcutter Shiny app: annotation codes, `prepare_results.R`, `run_leafviz.R`, the annotation-code mismatch |
| `references/jutils.md` | Tool-agnostic heatmaps, sashimi and Venn from rMATS/leafcutter/MntJULiP/MAJIQ output |
| `references/pygenometracks.md` | Multi-track figures: bedGraph/BigWig coverage, regtools junction arcs as BEDPE, tracks.ini |

## Reading Sashimi Plots (Interpretation Guide)

| Visual element | What it represents |
|----------------|--------------------|
| Filled coverage track | Read coverage at each genomic position (per sample, overlaid; the group mean with `-A mean`) |
| Arc / curve between exons | Junction-spanning reads; arc connects donor to acceptor |
| Number on arc | Junction-spanning alignment records (raw per sample; with `-A` the rounded group mean). ggsashimi skips only unmapped records: secondary and duplicate alignments are counted, so labels match a plain pysam CIGAR scan only when those are counted too (2 of 10 labels matched when secondary records were skipped) |
| Arc thickness | Often proportional to read count (tool-dependent) |
| Gene model below | Exons (boxes) and introns (lines) from GTF |
| Multiple parallel tracks | Per-sample (default) or per-group (with `-O`) |

**Junction count interpretation:** an arc's number is the count of alignments whose CIGAR has an `N` operation matching that intron. Higher = more usage. Compare inclusion vs skipping arcs to estimate PSI visually.

## Per-Tool Failure Modes

### ggsashimi: Off-Strand Junction Artifacts

**Trigger:** Stranded RNA-seq library plotted without strand specification.

**Mechanism:** ggsashimi reads BAM strand from CIGAR + flag; without strand info, antisense junctions appear as artifacts.

**Symptom:** Implausible junctions in regions with overlapping antisense genes; "noise" arcs at unexpected locations.

**Fix:** Set library strandedness with `-s MATE2_SENSE` (dUTP/TruSeq reverse-stranded PE; `-s MATE1_SENSE` for forward PE; verify orientation with RSeQC `infer_experiment.py`). `MATE1_SENSE`/`MATE2_SENSE` are paired-end only (`TypeError: ... 'NoneType' and 'bool'` on single-end); single-end uses `-s SENSE`/`ANTISENSE`. With `-s` the output is two files, `<prefix>_+.<fmt>` and `<prefix>_-.<fmt>`. Alternatively, pre-filter the BAM by strand with `samtools view -f 16` / `-F 16`.

### ggsashimi: Silent Failures (exit 0)

**Trigger:** BAM path typo in the TSV, region without reads, R package missing or ggplot2 error (e.g. `Unknown colour name`).

**Symptom:** rc 0 with a dropped sample, an empty figure, or no figure at all.

**Fix:** Check that every BAM exists, that the region has reads (`samtools view -c sample.bam chr1:100-200`), and that the figure file exists and is non-empty (recipes above).

## Best Practices

| Tip | Rationale |
|-----|-----------|
| Limit to 3-4 groups per figure | More becomes hard to read |
| Include 200-500 nt flanking exons | Show full splicing context |
| Check accessibility colors | Use ColorBrewer-safe palettes for color-blind readers |
| Always include a legend | Sashimi figures without legends are uninformative for non-experts |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| R error printed (`there is no package called ...`, `Unknown colour name`), **rc 0, no figure** | ggsashimi needs `R` on PATH with ggplot2/data.table/gridExtra/gtable (the samtools binary is not used; pysam is); or a bad palette colour | Install the R packages; assert the figure exists |
| `ggsashimi: ValueError: invalid contig` | Contig name differs between region and BAM (`chr1` vs `1`) | Map the name against the BAM header (recipes above) |
| `ERROR: No available bam files.` | Every path in the TSV is wrong (a single wrong path is dropped silently) | Check paths; relative paths resolve against the TSV's directory |
| `ERROR: Cannot apply aggregate function if overlay is not selected.` | `-A` without `-O` | Add `-O 3` |
| `RuntimeError: generator raised StopIteration` | `--shrink` with no junction passing `-M` | Lower `-M` or drop `--shrink` |
| `rmats2sashimiplot: unrecognized arguments: -t SE` | Old flag | `--event-type SE` |
| `Error: Must provide sample label and color for each entry in bam_files!`, rc 0, empty `Sashimi_plot/` | Colours per replicate given for 2 groups | Add `--group-info grouping.gf` or one colour per replicate |
| `pyGenomeTracks: InputError ... can not identify file type` | BAM track (unsupported) or missing `file_type` | Use bedGraph/BigWig; set `file_type` |
| `jutils.py venn-diagram: FileNotFoundError` on `a.tsv,b.tsv` | `--tsv-file-list` is a file of paths | Write the list file |
| `App dir must contain either app.R or server.R` | `run_leafviz.R` not started from the `leafviz/` directory | `cd leafcutter/leafviz` first |
| `prepare_results.R`: `<file> does not exist` | annotation_code prefix or a path is wrong | Re-run with correct paths |

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No junctions shown | `-M` above every junction's count (default is 1) | Lower `-M` |
| Plot too crowded | Many samples without aggregation | Use `-O 3` to overlay groups |
| Annotation missing or wrong gene | GTF lacks gene_name attribute or wrong build | Verify GTF version vs BAM reference; pre-filter the GTF to the relevant features |
| Memory issues on large regions | >100 kb regions with many samples | Plot smaller windows or pre-extract reads with samtools view |
| Y-axis dominated by one peak | Outlier sample | Filter the outlier out of the TSV |

## Related Skills

- differential-splicing - Identify events to plot; sashimi plots are validation
- splicing-quantification - Context for PSI values; sashimi provides visual confirmation
- data-visualization/genome-tracks - Multi-track figure design (pyGenomeTracks, Gviz)
- data-visualization/ggplot2-fundamentals - ggsashimi customization (extends ggplot2)
- data-visualization/color-palettes - Accessible color choices
- data-visualization/volcano-and-ma-plots - Volcano complement to sashimi
- data-visualization/heatmaps-clustering - Heatmap complement to sashimi

## References

- Katz et al 2010 *Nat Methods* - MISO sashimi plot original
- Garrido-Martin et al 2018 *PLoS Comput Biol* - ggsashimi
- Yang et al 2021 *Bioinformatics* - Jutils
- Vaquero-Garcia et al 2016 *eLife* - MAJIQ / VOILA
- Li et al 2018 *Nat Genet* - leafcutter / leafviz
- Ramirez et al 2018 *Nat Commun* - pyGenomeTracks
