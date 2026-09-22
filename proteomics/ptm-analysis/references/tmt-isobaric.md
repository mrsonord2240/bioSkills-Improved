# TMT / Isobaric Plexes (MSstatsPTM)

### TMT / isobaric plexes -- same adjustment, three different calls

**Goal:** The same protein-adjusted site testing when the enriched and global runs are isobaric-labelled plexes (the dominant platform for large-cohort phosphoproteomics).

**Approach:** Only three things change. `labeling_type = 'TMT'` on the converter (default `'LF'`), `dataSummarizationPTM_TMT` instead of `dataSummarizationPTM` -- a different FUNCTION, not a flag -- and `data.type = 'TMT'` in `groupComparisonPTM`. The annotation is the MSstatsTMT one: `Run`, `Fraction`, `TechRepMixture`, `Channel`, `Condition`, `Mixture`, `BioReplicate` (`Raw.file` is also accepted). Mismatching the two halves fails loudly in BOTH directions, but neither message says "wrong labeling type": a label-free annotation with `labeling_type = 'TMT'` gives *`Extra columns included in the annotation file that are not required ... Run, Raw.file, Fraction, TechRepMixture, Channel, Condition, Mixture, BioReplicate`* -- read that list as the spec for the TMT annotation -- and TMT evidence left on the default `'LF'` gives one of several equally unhelpful errors (see Common Errors). EVERY plex needs a pooled reference channel carried as `Condition = 'Norm'`; `reference_norm = TRUE` (default) uses it to put plexes on a common scale -- the IRS bridge -- and `remove_norm_channel = TRUE` (default) drops it before testing. Without one, cross-plex comparison is invalid and no amount of downstream modelling repairs it.

Two TMT-specific traps for the ADJUSTMENT itself, not just for the quant:

1. **Ratio compression biases the subtraction, it does not merely attenuate it.** Co-isolated precursors pull every reporter ratio toward 1, and the PTM and PROTEIN datasets are not compressed equally -- enriched phosphopeptide runs are a sparser, differently-interfered precursor space than the global run. So `dFC_adj = dFC_PTM - dFC_protein` subtracts two differently-shrunken numbers. Use MS3/SPS or FAIMS to reduce it, and read TMT effect sizes as lower bounds; the SIGN and ranking survive, the magnitude does not.
2. **Label the enriched and global aliquots in the SAME plex** where the design allows. Split across plexes, each dataset carries its own reference-channel normalization and the site-to-protein subtraction inherits both.

```bash
Rscript scripts/msstatsptm_tmt.R dir=<data_dir> out=<out_dir>   # writes <out_dir>/adjusted_sites_tmt.csv
```

`scripts/msstatsptm_tmt.R` carries the whole route: the class-I pre-filter on the enriched evidence (as in the label-free script), `MaxQtoMSstatsPTMFormat(labeling_type = 'TMT', ...)` with the global run, `dataSummarizationPTM_TMT`, an explicit `Treatment vs Control` contrast, `groupComparisonPTM(data.type = 'TMT')`, and the site rows of `ADJUSTED.Model`. Its comments hold the `Channel` 0-indexing rule (annotation `channel.0` .. `channel.9` for a 10-plex, read off the evidence header) and the `Condition = 'Norm'` reference channel. From the adjusted table on, the TREAT-style threshold, the PTM.Model-vs-ADJUSTED.Model comparison and the KSEA step are identical to the label-free route.

Starting from `Phospho (STY)Sites.txt` instead of `evidence.txt` (the `sites_data =` argument) is the one place `TMT_keyword` matters: there the converter builds column names as `Reporter.intensity.corrected.<n>.<TMT_keyword><mixture>`, so `TMT_keyword` must match how the site table's reporter columns were named. It is ignored on the `evidence =` route shown above.
