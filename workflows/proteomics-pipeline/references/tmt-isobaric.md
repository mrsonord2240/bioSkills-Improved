# TMT / iTRAQ isobaric labeling

Read this when the data are TMT/iTRAQ: reporter extraction from mzML, CoA impurity correction with the orientation guard, and the multi-plex MSstatsTMT route with the reference-channel bridge.

### TMT/iTRAQ Isobaric Labeling
Reporter extraction is a spectra-level step, not a text-matrix read. Within a single plex the channels are co-isolated/co-fragmented in the same MS2 event, so relative ratios are stable; but MULTI-batch TMT CANNOT be compared across plexes without an IRS bridge (a pooled reference channel in every plex; Plubell 2017). Route to proteomics/quantification for the mechanics.
```r
library(MSnbase)

# Extract reporter ions from spectra (NOT readMSnSet, which loads an existing text matrix)
raw <- readMSData('tmt.mzML', mode = 'onDisk')
tmt_data <- quantify(raw, reporters = TMT10, method = 'max')
# Correct isobaric impurity cross-talk with the LOT-SPECIFIC matrix from the reagent CoA.
# makeImpuritiesMatrix(x = 10, edit = FALSE) returns a MANUFACTURER TEMPLATE, not an identity
# matrix -- its diagonal runs 0.928-0.965 and 5% of 126 lands in 127C. It is a shape check only;
# the numbers are lot-specific. (edit = TRUE, the default, opens an editor and blocks in scripts.)
# Do NOT use makeImpuritiesMatrix(filename = ...) for TMT10/TMTpro: it reads a CoA laid out by
# Da OFFSET and places each column k POSITIONS away in the reporter list, which is only correct
# for non-interleaved reagents (TMT6, iTRAQ). TMT10/TMTpro interleave N and C, so the +1 Da
# neighbour of 126 is 127C -- TWO positions away -- and the filename route silently writes the
# bleed into 127N instead. Build the matrix by CHANNEL NAME and hand it to purityCorrect:
# tmt10_coa.csv: a square percentage matrix, rows = SOURCE reagent, columns = OBSERVED channel,
# both labelled with the channel names (126, 127N, 127C, ...); diagonal = the lot's purity.
coa <- as.matrix(read.csv('tmt10_coa.csv', row.names = 1, check.names = FALSE))
stopifnot(nrow(coa) == ncol(coa), setequal(rownames(coa), reporterNames(TMT10)))
coa <- coa[reporterNames(TMT10), reporterNames(TMT10)]   # force the quant's channel order
# ORIENTATION CHECK. A transposed sheet has the right channel names, the right shape and produces
# NO negative values, so the negative-count check below never sees it -- yet it makes the
# correction 2.7x worse than the correct orientation (median relative error 0.0015 vs 0.0006 on
# the TMT10 template), still better than doing nothing and therefore silent. Test the orientation
# directly: a ROW is one reagent's isotopic envelope and sums to 100% by construction (minus what
# falls off the ends of the channel list); a COLUMN sums over different reagents and has no such
# constraint. If the columns fit 100 better than the rows, the sheet is the wrong way round.
row_dev <- sum((rowSums(coa) - 100)^2); col_dev <- sum((colSums(coa) - 100)^2)
if (col_dev < row_dev)
    stop('CoA looks TRANSPOSED: column sums fit 100% better than row sums (', round(col_dev, 1),
         ' vs ', round(row_dev, 1), '). Rows must be the SOURCE reagent, columns the OBSERVED ',
         'channel -- transpose the sheet or re-read the lot certificate.')
stopifnot(all(diag(coa) == apply(coa, 1, max)),   # each reagent's own channel must dominate its row
          all(diag(coa) > 50))                    # a CoA is percentages; < 50 means fractions were read
impurities <- coa / 100                                  # CoA percentages -> fractions
tmt_data <- purityCorrect(tmt_data, impurities)
stopifnot(sum(exprs(tmt_data) < 0, na.rm = TRUE) == 0)   # negatives = a grossly wrong matrix (NOT a transposition test; see above)

```

Multi-plex TMT (the common case: two or more plexes) -- do NOT concatenate plexes directly. MSstatsTMT applies the reference-channel (IRS-style) bridge during summarization. MaxQuant route, checked on MSstatsTMT 2.14.2 against its bundled 5-plex `evidence` / `proteinGroups` / `annotation.mq`.
```r
library(MSstatsTMT)

# Multi-plex TMT from MaxQuant. Each plex (Mixture) carries a pooled reference channel, annotated
# Condition = 'Norm' (MSstatsTMT requires that exact label); that channel is the bridge. annotation.csv has one row per (Run, Channel) with
# columns Run, Fraction, TechRepMixture, Channel, Condition, Mixture, BioReplicate.
# quote = '' / comment.char = '' as in references/msstats.md.
evidence <- read.table('evidence.txt', sep = '\t', header = TRUE, quote = '', comment.char = '')
proteinGroups <- read.table('proteinGroups.txt', sep = '\t', header = TRUE, quote = '', comment.char = '')
stopifnot(nrow(evidence) == length(readLines('evidence.txt')) - 1,
          nrow(proteinGroups) == length(readLines('proteinGroups.txt')) - 1)
annotation <- read.csv('annotation.csv')
stopifnot('Norm' %in% annotation$Condition)   # no reference channel = nothing to bridge plexes with

tmt_input <- MaxQtoMSstatsTMTFormat(evidence, proteinGroups, annotation, use_log_file = FALSE,
                                     verbose = FALSE)

# Summarize to protein level. Within-plex global median normalization is followed by reference-channel
# normalization: every plex is rescaled to its own 'Norm' channel, which is the cross-plex (IRS-style)
# bridge, and the Norm channel is then dropped. Do not concatenate plexes before this step.
summ <- proteinSummarization(tmt_input, method = 'msstats', global_norm = TRUE, reference_norm = TRUE,
                             remove_norm_channel = TRUE, use_log_file = FALSE, verbose = FALSE)

# Contrasts: one row per comparison, columns = the Condition levels that survive (Norm is removed),
# in sorted order -- built from the levels, against the first one, as in the limma block in SKILL.md.
lv <- sort(setdiff(unique(as.character(annotation$Condition)), 'Norm'))
comparison <- t(sapply(lv[-1], function(l) as.numeric(lv == l) - as.numeric(lv == lv[1])))
colnames(comparison) <- lv
rownames(comparison) <- paste0(lv[-1], '_vs_', lv[1])
# moderated = TRUE borrows variance across proteins (limma-style); groupComparisonTMT adjusts within
# each contrast, so adjust ACROSS the rows yourself when you report several (as in references/msstats.md).
tmt_res <- groupComparisonTMT(summ, contrast.matrix = comparison, moderated = TRUE,
                              adj.method = 'BH', use_log_file = FALSE, verbose = FALSE)$ComparisonResult
tmt_res$adj.pvalue.global <- p.adjust(tmt_res$pvalue, method = 'BH')
print(table(tmt_res$Label, tmt_res$adj.pvalue < 0.05))
```
