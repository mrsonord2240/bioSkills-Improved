## pyGenomeTracks for Multi-Track Figures

**Goal:** Combine splicing with chromatin or coverage tracks for publication figures.

**Approach:** Build coverage bedGraphs and a junction BEDPE from the BAMs, define tracks in an INI file (genes, bedGraph/BigWig, BED, links), then run `pyGenomeTracks --tracks tracks.ini --region ... -o figure.pdf`. pyGenomeTracks 3.9 cannot draw a BAM (`InputError ... can not identify file type`).

```bash
# one merged BAM per group; -split is essential: without it introns are filled with coverage
samtools merge -f ctrl_merged.bam ctrl1.bam ctrl2.bam ctrl3.bam && samtools index ctrl_merged.bam
samtools merge -f trt_merged.bam trt1.bam trt2.bam trt3.bam && samtools index trt_merged.bam
bedtools genomecov -ibam ctrl_merged.bam -split -bga > ctrl.bedgraph
bedtools genomecov -ibam trt_merged.bam -split -bga > trt.bedgraph
```

```ini
[gene_models]
file = annotation.gtf
height = 3
title = GENCODE v45
fontsize = 10
file_type = gtf

[ctrl_coverage]
file = ctrl.bedgraph
title = Control
color = #1f77b4
height = 3
min_value = 0
max_value = 200
file_type = bedgraph

[trt_coverage]
file = trt.bedgraph
title = Treatment
color = #ff7f0e
height = 3
min_value = 0
max_value = 200
file_type = bedgraph

[junctions]
file = junctions.bedpe
title = Junctions
height = 5
file_type = links
links_type = arcs
```

Arc height grows with the junction's span, so a short `[junctions]` track crops the widest arc (`height = 2` did on a 3-exon locus whose skipping junction spans 60% of the window; 5 shows it whole): raise `height`, or narrow `--region`, until the widest arc is complete, and look at the figure.

Tracks are scaled independently: set the same `min_value`/`max_value` on both coverage tracks (pick `max_value` from the data) or the two groups are not comparable. A BigWig made from the same `-split` bedGraph works too (`file_type = bigwig`).

The `junctions.bedpe` file must be in **BEDPE format** (6 columns: chr1 start1 end1 chr2 start2 end2 [+ optional score]). Convert from regtools .bed12 junctions (the score is the read count summed over the merged BAM; `-s XS` needs XS-tagged BAMs, otherwise the strand is `?`):

```bash
samtools merge -f all_merged.bam ctrl_merged.bam trt_merged.bam && samtools index all_merged.bam   # regtools needs an indexed BAM
regtools junctions extract -s XS -o regtools_junctions.bed all_merged.bam
# regtools BED12 column 11 is blockSizes (anchor_left, anchor_right);
# column 12 is blockStarts (0, intron_length + anchor_left).
# Intron start = chromStart + anchor_left = $2 + a[1]
# Intron end   = chromStart + blockStarts[2] = $2 + b[2]
awk 'BEGIN{OFS="\t"} {split($11,a,","); split($12,b,","); s=$2+a[1]; e=$2+b[2]; print $1, s, s+1, $1, e-1, e, $5}' \
    regtools_junctions.bed > junctions.bedpe
```

```bash
pyGenomeTracks --tracks tracks.ini --region chr17:43094000-43125000 -o figure.pdf
```
