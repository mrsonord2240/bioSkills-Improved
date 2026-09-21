## rmats2sashimiplot

**Goal:** Plot directly from rMATS event coordinates without manual region calculation.

**Approach:** Filter the rMATS event file to the events to plot (it draws every row), pass BAM lists + a group file + event type, then check the output: rmats2sashimiplot **exits 0 when it fails** and leaves `Sashimi_plot/` empty. The contig is matched to the BAM header automatically (`chrX` in the event file, `X` in the BAM works).

```bash
# rmats2sashimiplot plots every row: keep only significant events (columns found by header name)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next}
    $c["FDR"]<0.05 && ($c["IncLevelDifference"]>0.1 || $c["IncLevelDifference"]<-0.1)' \
    rmats_output/SE.MATS.JC.txt > sig.SE.MATS.JC.txt

# group file: "label: first-last", 1-based over the --b1 replicates then the --b2 replicates
printf 'Control: 1-3\nTreatment: 4-6\n' > grouping.gf

rmats2sashimiplot \
    --b1 ctrl1.bam,ctrl2.bam,ctrl3.bam \
    --b2 trt1.bam,trt2.bam,trt3.bam \
    --event-type SE \
    -e sig.SE.MATS.JC.txt \
    --l1 Control \
    --l2 Treatment \
    -o sashimi_rmats \
    --exon_s 1 \
    --intron_s 5 \
    --group-info grouping.gf \
    --color '#1f77b4,#ff7f0e'

n_events=$(( $(wc -l < sig.SE.MATS.JC.txt) - 1 ))
n_pdf=$(find sashimi_rmats/Sashimi_plot -name '*.pdf' -size +0 2>/dev/null | wc -l)
[ "$n_pdf" -eq "$n_events" ] || { echo "rmats2sashimiplot wrote $n_pdf of $n_events figures" >&2; exit 1; }
```

`--event-type` (4.0.0; the old `-t SE` is rejected, rc 2) takes SE, A5SS, A3SS, MXE or RI. `--exon_s 1 --intron_s 5` draws introns at 1/5 of their real length. `--group-info` gives one plot per group (arc labels = group mean, plus the group's mean IncLevel); without it there is one plot per replicate and `--color` needs one colour per replicate, otherwise it prints `Error: Must provide sample label and color for each entry in bam_files!` and still exits 0.
