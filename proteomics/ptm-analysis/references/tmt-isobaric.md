# TMT / Isobaric Plexes (MSstatsPTM)

### TMT / isobaric plexes -- same adjustment, three different calls

**Goal:** The same protein-adjusted site testing when the enriched and global runs are isobaric-labelled plexes (the dominant platform for large-cohort phosphoproteomics).

**Approach:** Only three things change. `labeling_type = 'TMT'` on the converter (default `'LF'`), `dataSummarizationPTM_TMT` instead of `dataSummarizationPTM` -- a different FUNCTION, not a flag -- and `data.type = 'TMT'` in `groupComparisonPTM`. The annotation is the MSstatsTMT one: `Run`, `Fraction`, `TechRepMixture`, `Channel`, `Condition`, `Mixture`, `BioReplicate` (`Raw.file` is also accepted). Mismatching the two halves fails loudly in BOTH directions, but neither message says "wrong labeling type": a label-free annotation with `labeling_type = 'TMT'` gives *`Extra columns included in the annotation file that are not required ... Run, Raw.file, Fraction, TechRepMixture, Channel, Condition, Mixture, BioReplicate`* -- read that list as the spec for the TMT annotation -- and TMT evidence left on the default `'LF'` gives one of several equally unhelpful errors (see Common Errors). EVERY plex needs a pooled reference channel carried as `Condition = 'Norm'`; `reference_norm = TRUE` (default) uses it to put plexes on a common scale -- the IRS bridge -- and `remove_norm_channel = TRUE` (default) drops it before testing. Without one, cross-plex comparison is invalid and no amount of downstream modelling repairs it.

Two TMT-specific traps for the ADJUSTMENT itself, not just for the quant:

1. **Ratio compression biases the subtraction, it does not merely attenuate it.** Co-isolated precursors pull every reporter ratio toward 1, and the PTM and PROTEIN datasets are not compressed equally -- enriched phosphopeptide runs are a sparser, differently-interfered precursor space than the global run. So `dFC_adj = dFC_PTM - dFC_protein` subtracts two differently-shrunken numbers. Use MS3/SPS or FAIMS to reduce it, and read TMT effect sizes as lower bounds; the SIGN and ranking survive, the magnitude does not.
2. **Label the enriched and global aliquots in the SAME plex** where the design allows. Split across plexes, each dataset carries its own reference-channel normalization and the site-to-protein subtraction inherits both.

```r
library(MSstatsPTM)
rd <- function(f) read.table(f, sep = '\t', header = TRUE, quote = '')

# Class-I pre-filter on the ENRICHED evidence, exactly as in the label-free route above.
ev <- rd('evidence_phospho_tmt.txt')
site_prob <- vapply(regmatches(ev$Phospho..STY..Probabilities,
                               gregexpr('(?<=\\()[0-9.]+(?=\\))', ev$Phospho..STY..Probabilities, perl = TRUE)),
                    function(p) if (length(p)) max(as.numeric(p)) else NA_real_, numeric(1))
ev <- ev[grepl('Phospho \\(STY\\)', ev$Modified.sequence) & !is.na(site_prob) & site_prob >= 0.75, ]

# TMT annotation: Run, Fraction, TechRepMixture, Channel, Condition, Mixture, BioReplicate.
# Channel names follow the reporter-column suffixes MaxQuant wrote, and those are 0-indexed: a
# 10-plex has 'Reporter intensity corrected 0' .. '9', so the annotation needs 'channel.0' .. 'channel.9'
# (the CORRECTED columns, not the raw reporters). Read the suffixes off your own evidence header;
# 'channel.1' .. 'channel.10' is rejected with 'the channel name must be matched with that in input
# data', which never mentions the off-by-one.
# Give the pooled reference channel Condition = 'Norm' in EVERY plex.
input <- MaxQtoMSstatsPTMFormat(
  evidence = ev,
  annotation = read.csv('annotation_ptm_tmt.csv'),
  fasta_path = 'uniprot_human.fasta',
  evidence_prot = rd('evidence_global_tmt.txt'),
  proteinGroups = rd('proteinGroups_global.txt'),
  annotation_protein = read.csv('annotation_protein_tmt.csv'),
  labeling_type = 'TMT',          # the single converter switch; default is 'LF'
  mod_id = '\\(Phospho \\(STY\\)\\)',
  which_proteinid_ptm = 'Proteins',
  which_proteinid_protein = 'Proteins',
  use_unmod_peptides = FALSE)
stopifnot('PROTEIN' %in% names(input))

# TMT summarization is its own function. reference_norm / reference_norm.PTM (default TRUE) apply
# the 'Norm'-channel bridge; remove_norm_channel (default TRUE) drops that channel afterwards, so
# the contrast below names only the biological conditions.
summarized <- dataSummarizationPTM_TMT(input, use_log_file = FALSE, append = FALSE)

contrast <- matrix(c(-1, 1), nrow = 1, dimnames = list('Treatment vs Control', c('Control', 'Treatment')))
result <- groupComparisonPTM(summarized, data.type = 'TMT', contrast.matrix = contrast)

adjusted <- result$ADJUSTED.Model
adjusted <- adjusted[grepl('_[STY][0-9]+', adjusted$Protein), ]
# From here the TREAT-style threshold, the PTM.Model-vs-ADJUSTED.Model comparison and the KSEA
# block below are identical to the label-free route -- only the three calls above differ.
```

Starting from `Phospho (STY)Sites.txt` instead of `evidence.txt` (the `sites_data =` argument) is the one place `TMT_keyword` matters: there the converter builds column names as `Reporter.intensity.corrected.<n>.<TMT_keyword><mixture>`, so `TMT_keyword` must match how the site table's reporter columns were named. It is ignored on the `evidence =` route shown above.
