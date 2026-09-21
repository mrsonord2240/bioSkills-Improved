## leafviz Shiny App

**Goal:** Browse leafcutter clusters with intron-level effects and sashimi-like plots.

**Approach:** leafviz is a script directory inside the leafcutter repo (not an R package: `library(leafviz)` and `run_leafviz()` do not exist). Build annotation files from the GTF, prepare the results `.RData`, then launch the Shiny app from the `leafviz` directory. leafcutter is a GitHub R package (`devtools::install_github('davidaknowles/leafcutter/leafcutter')`), not Bioconductor (the as-shipped 0.2.9 fails to build against rstan >= 2.33 because of the old Stan array syntax; confirm `library(leafcutter)` loads); the `leafviz/` directory comes with its repo. Checked end to end on the planted 3v3 leafcutter results (leafcutter 0.2.9; the app serves HTTP 200).

```bash
# annotation_code = prefix of four files (_all_exons.txt.gz, _all_introns.bed.gz, _fiveprime.bed.gz, _threeprime.bed.gz);
# build them from the GTF version used in the differential analysis
perl leafcutter/leafviz/gtf2leafcutter.pl -o annot annotation.gtf

# groups.txt = the support file given to leafcutter_ds.R (sample <TAB> condition)
Rscript leafcutter/leafviz/prepare_results.R \
    -o leafviz.RData \
    -m groups.txt \
    leafcutter_perind_numers.counts.gz \
    ds_results_cluster_significance.txt \
    ds_results_effect_sizes.txt \
    annot

# runApp() uses the working directory: start from leafviz/, pass the .RData by absolute path
cd leafcutter/leafviz && Rscript run_leafviz.R /abs/path/leafviz.RData    # prints "Listening on http://127.0.0.1:<port>"
```

`download_human_annotation_codes.sh` in the same directory fetches prebuilt hg19 codes. Useful for cohort-level interactive filtering of clusters.

### leafviz: Annotation Codes Mismatch

**Trigger:** Using leafviz with annotation_codes from different GENCODE version than leafcutter clusters.

**Mechanism:** annotation_codes encodes intron-to-event-class mapping per GTF version.

**Symptom:** Many clusters show as "unannotated" despite being in canonical GTF.

**Fix:** Generate annotation_codes with `gtf2leafcutter.pl` from the same GTF used in differential analysis.
