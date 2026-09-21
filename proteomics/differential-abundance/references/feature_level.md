# Feature-Level Testing (msqrob2, MSstats) -- Peptide Table In, Protein Calls Out

Loaded from `SKILL.md` when a peptide/precursor table exists. Checked 2026-09-21 on msqrob2 1.14.1, QFeatures 1.16.0, MsCoreUtils 1.18.0, MSstats 4.14.2 (R 4.4.3 / Bioconductor 3.20). Every block chains through the ambient objects `sample_info`, `peptide_wide`, `evidence`, `protein_groups` and `annotation`; the end-to-end version is `examples/msqrob2_peptide_level.R`.

## msqrob2 Workflow (R)

**Goal:** Test from the peptide/precursor table itself (MaxQuant `evidence.txt`, a DIA-NN report, a PSM table) so the between-peptide spread and the number of observations set the standard error, instead of collapsing to one number per protein per run first.

**Approach:** Build a `QFeatures` object from a wide peptide matrix, log-transform, keep peptides seen in at least 2 runs, aggregate to protein with `robustSummary` (Huber M-estimation, which is what downweights an outlier peptide), then fit `msqrob` -- ridge/robust regression with empirical-Bayes variance moderation -- and test the contrast. `hypothesisTest` writes its result into `rowData(pe[['protein']])`, one data frame per contrast. Two things the object does silently and you must undo: proteins with no observation in a condition come back with `adjPval = NA` rather than in a list, and any per-run normalization you apply here is applied to a *different peptide set in each run* (see the centring section). Report the untestable and undetected proteins separately, as `SKILL.md` requires. Runnable end to end: `examples/msqrob2_peptide_level.R` (takes `evidence.txt` plus an annotation file, or simulates a peptide table when given no arguments).

```r
library(QFeatures)
library(msqrob2)

# peptide_wide: one row per precursor; columns 'feature', 'protein', then one intensity column per run
runs <- sample_info$run
col_data <- data.frame(quantCols = runs, condition = factor(sample_info$condition),
                       sample = factor(runs), row.names = runs)  # quantCols column is required by readQFeatures
pe <- readQFeatures(assayData = peptide_wide, quantCols = runs, colData = col_data, name = 'peptideRaw')
pe <- zeroIsNA(pe, 'peptideRaw')
pe <- logTransform(pe, base = 2, i = 'peptideRaw', name = 'peptideLog')
rowData(pe[['peptideLog']])$nNonZero <- rowSums(!is.na(assay(pe[['peptideLog']])))
pe <- filterFeatures(pe, ~ nNonZero >= 2, keep = TRUE)  # keep=TRUE: the variable exists only on peptideLog

# undetected list, taken BEFORE aggregation: msqrob2 reports these as adjPval = NA, not as a list
cond <- colData(pe)$condition  # colData lives on the QFeatures object; pe[['peptideLog']]$condition is NULL
obs <- sapply(levels(cond), function(g)
  tapply(rowSums(!is.na(assay(pe[['peptideLog']])[, cond == g, drop = FALSE])),
         rowData(pe[['peptideLog']])$protein, sum))
undetected <- rownames(obs)[apply(obs, 1, min) == 0]  # report as "undetected in group X", never as a fold change

pe <- aggregateFeatures(pe, i = 'peptideLog', fcol = 'protein', name = 'protein',
                        fun = MsCoreUtils::robustSummary, na.rm = TRUE)
pe <- msqrob(pe, i = 'protein', formula = ~condition, robust = TRUE)
L <- makeContrast('conditionTreatment = 0', parameterNames = 'conditionTreatment')
pe <- hypothesisTest(pe, i = 'protein', contrast = L)

res <- rowData(pe[['protein']])$conditionTreatment  # columns: logFC, se, df, t, pval, adjPval (no adj.P.Val)
res$protein <- rownames(pe[['protein']])
untestable <- res$protein[is.na(res$adjPval)]       # report these; they are not "not significant"
res <- res[!is.na(res$adjPval), ]
```

To keep every peptide as its own degree of freedom instead of summarizing first, fit the mixed model over the peptide assay. `sample` must be a column of `colData` and `feature` a column of `rowData`:

```r
pe <- msqrobAggregate(pe, i = 'peptideLog', fcol = 'protein', name = 'proteinLmer',
                      formula = ~condition + (1 | sample) + (1 | feature), ridge = FALSE)
pe <- hypothesisTest(pe, i = 'proteinLmer', contrast = L)
```

`ridge = TRUE` needs MORE than two parameters in the mean model, so it is refused outright on a plain two-group comparison ("The mean model must have more than two parameters for ridge regression"); use `ridge = FALSE`, or drop the intercept with `~ -1 + condition` as the error message itself suggests. Ridge is for three or more groups, or a group factor plus covariates.

## MSstats Workflow (R)

**Goal:** Test a feature-level design with run/subject structure -- technical replicates, nested or repeated measures, SRM/PRM/DIA -- from the search engine's own output tables.

**Approach:** Convert with the vendor-specific importer, summarize with `dataProcess`, then test the contrast with `groupComparison`. The contrast matrix is columns-by-condition and its `dimnames` must match `levels(proc$ProteinLevelData$GROUP)` exactly. Keep `MBimpute = FALSE` unless you want the AFT censored imputation: it is the only place MSstats imputes, and the Skill's position is to model dropout rather than fill it. On a 4 v 4 label-free set the two settings differed by 0.4 percentage points of realized FDR, so the AFT imputation is not buying accuracy either. Proteins present in only one condition come back with an infinite `log2FC` and `issue == 'oneConditionMissing'`; split them out and report them as undetected.

**The `normalization` argument is the per-run median normalization the centring section below warns about.** `'equalizeMedians'` is MSstats's default. On the 4 v 4 audit peptide set it gave median log2FC -0.184, median null-protein SE 0.119 and 101 calls with 21 false positives (20.8% realized FDR at nominal 5%); `normalization = FALSE` gave median log2FC +0.008, SE 0.204 and 79 calls with 0 false positives. Keep `'equalizeMedians'` only if the centring checks below pass; otherwise pass `FALSE` and normalize at the protein level after summarization (proteomics/quantification).

```r
library(MSstats)

input <- MaxQtoMSstatsFormat(evidence = evidence, proteinGroups = protein_groups,
                             annotation = annotation, use_log_file = FALSE)  # annotation: Raw.file, Condition, BioReplicate, IsotopeLabelType
proc <- dataProcess(input,
                    normalization = 'equalizeMedians',  # per-run median normalization: run the centring checks; FALSE = none (audit set: -0.184 / 20.8% FDR vs +0.008 / 0.0%)
                    summaryMethod = 'TMP', censoredInt = 'NA', MBimpute = FALSE, use_log_file = FALSE)

contrast <- matrix(c(-1, 1), nrow = 1,
                   dimnames = list('Treatment-Control', c('Control', 'Treatment')))  # order = levels(GROUP)
res <- groupComparison(contrast.matrix = contrast, data = proc, use_log_file = FALSE)$ComparisonResult
res$Protein <- as.character(res$Protein)

undetected <- res[res$issue %in% 'oneConditionMissing', ]   # infinite log2FC; report as undetected, not as a ratio
tested <- res[is.finite(res$log2FC) & !is.na(res$adj.pvalue), ]
# columns: Protein, Label, log2FC, SE, Tvalue, DF, pvalue, adj.pvalue, issue (adj.pvalue is the BH p)
```

## Centring checks -- run both before reading any feature-level result

**Goal:** Catch the failure that turns a feature-level test anticonservative across the whole table.

**Approach:** Per-run median normalization of a peptide table does two things at once, and a feature-level test is exposed to both. (1) Each run's median is taken over the peptides detected in THAT run, and the detected sets differ, so equalizing the medians transfers a detection-composition difference into a uniform between-condition offset. (2) It also removes real run-to-run loading variation within each condition, which shrinks the residual the test divides by. Neither is enough alone, measured on one 4 v 4 label-free peptide set with planted truth (msqrob2, nominal FDR 5%):

| variant | median log2FC | null-protein median SE | calls | false positives |
|---|---|---|---|---|
| no normalization | -0.005 | 0.206 | 80 | 0 |
| uniform -0.20 injected into every treatment run (offset only) | -0.205 | 0.206 | 79 | 0 |
| each run centred on its own condition's median (SE shrink only) | -0.016 | 0.116 | 81 | 1 |
| per-run median over all detected peptides (both) | -0.199 | 0.116 | 111 | 31 (27.9%) |

The offset is the detectable symptom of the pair, not the cause on its own, so there are two cheap checks: the offset of the result table, and the residual SD of the peptide table before and after the normalization. The residual SD fell to 0.76 of its un-normalized value under per-run median centring (both the all-peptide and the within-condition-only variants) and stayed at 1.00 under the uniform offset. Fix a failure upstream (proteomics/quantification: normalize on features present in every run, or at the protein level after summarization), not by re-centring the p-values.

```r
# Trip-wires, not distributional bounds: both constants were read off ONE 4 v 4 label-free set.
OFFSET_MAX   <- 0.05  # |median log2FC| over the tested proteins
SD_RATIO_MIN <- 0.85  # residual SD after / before per-run normalization (0.76 where it produced 31 false positives)

# median within-condition SD over peptides (log2 matrix, features x runs; NA-tolerant)
resid_sd <- function(L, cond) {
  ss <- df <- 0
  for (g in levels(cond)) {
    m <- L[, cond == g, drop = FALSE]
    ss <- ss + rowSums((m - rowMeans(m, na.rm = TRUE))^2, na.rm = TRUE)
    df <- df + pmax(rowSums(!is.na(m)) - 1, 0)
  }
  median(sqrt(ss[df > 0] / df[df > 0]))
}

# 1. Residual check, only if you normalize. msqrob2 route: run this after filterFeatures and before
#    aggregateFeatures, then aggregate and test from 'peptideNorm' (aggregateFeatures(pe, i = 'peptideNorm', ...)).
#    MSstats route: fit twice, normalization = FALSE and the default, and compare median(tested$SE)
#    instead: same warning, same SD_RATIO_MIN.
pe <- normalize(pe, i = 'peptideLog', name = 'peptideNorm', method = 'center.median')
sd_ratio <- resid_sd(assay(pe[['peptideNorm']]), cond) / resid_sd(assay(pe[['peptideLog']]), cond)
if (sd_ratio < SD_RATIO_MIN) warning(sprintf(
  'residual SD fell to %.2f of its un-normalized value: the normalization is removing within-condition variance and shrinking the SEs. Check the offset below and compare with normalization off.',
  sd_ratio))

# 2. Offset check on the result table (msqrob2: res$logFC; MSstats: tested$log2FC)
offset <- median(res$logFC, na.rm = TRUE)
if (abs(offset) > OFFSET_MAX) stop(sprintf(
  'median log2FC = %+.3f: the contrast is not centred. Re-normalize (proteomics/quantification) before reading this table.',
  offset))
```

Passing both checks does not certify the normalization; it only means these two symptoms are absent. The checks are also imperfect in both directions on the audit set: the offset stop would have blocked a harmless uniform -0.20 shift, and the within-condition-only variant passes the offset check (-0.016) and is caught only by the SD ratio. A drop in residual SD is also what normalization is FOR when the loading variation is real, so read a warning together with the offset. A study that genuinely expects a global shift must override `OFFSET_MAX` deliberately.

## Failure mode: per-run median normalization under a feature-level test

**Trigger:** `dataProcess(normalization = 'equalizeMedians')`, `normalize(method = 'center.median')` or any per-run median centring applied to the PEPTIDE table, then msqrob2/MSstats.
**Mechanism:** the offset transfer plus the residual shrinkage described above. A protein-summary test would mostly shrug off the offset; a feature-level test has standard errors small enough to call it significant once the residual has also been shrunk.
**Symptom:** `median(log2FC)` over all tested proteins far from 0 while the raw peptide table's is ~0; the extra calls are low-|logFC|, low-SE proteins; the true hits are recovered either way (48 up / 32 down in every msqrob2 variant).
**Fix:** run the centring checks; normalize on features present in every run, or after summarization at the protein level. Do not re-centre p-values.

## Threshold

| Threshold | Source | Rationale |
|---|---|---|
| `OFFSET_MAX` 0.05, `SD_RATIO_MIN` 0.85 | measured, one 4 v 4 set | realized FDR at nominal 5% was 0.0% at offset -0.005, 8.0% at -0.107 (median over complete-case peptides), 27.9% at -0.199; empirical trip-wires, not distributional bounds |
| msqrob2 ridge needs > 2 mean-model parameters | msqrob2 1.14.1 | a two-group `~condition` mean model has 2, so `ridge = TRUE` errors; `ridge = FALSE` or `~ -1 + condition` |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| a block of low-\|logFC\| calls from msqrob2/MSstats and `median(log2FC)` far from 0 | per-run median normalization on the peptide table | see the centring checks |
| `readQFeatures`: `'colData' must contain a column called 'quantCols'` | QFeatures 1.16 requires the run names in `colData$quantCols` | `data.frame(quantCols = runs, condition = ..., row.names = runs)` |
| msqrob2: `The mean model must have more than two parameters for ridge regression` | `ridge = TRUE` on a two-group design | `ridge = FALSE`, or drop the intercept (`~ -1 + condition`) |
| msqrob2: `Variable sample is not found in coldata or rowdata` | `msqrobAggregate` formula names a term that is in neither | add `sample` to `colData`; `feature` must be a `rowData` column of the peptide assay |
| `pe[['peptideLog']]$condition` is `NULL` | in QFeatures 1.16 `colData` is held on the QFeatures object, not on each assay | `colData(pe)$condition` (or `getWithColData(pe, i)` for a standalone assay) |
| msqrob2 results have no `adj.P.Val`; some rows are all `NA` | `hypothesisTest` writes `logFC, se, df, t, pval, adjPval` into `rowData`, and returns `NA` for proteins it cannot fit | read `adjPval`; report the `is.na(adjPval)` rows separately, not as non-significant |
| MSstats rows with infinite `log2FC` | `issue == 'oneConditionMissing'` | split them out and report as undetected; never as a ratio |
| MaxQuant `evidence.txt`: the `Reverse`/`Potential contaminant` filter drops every row | empty flag columns are read as logical `NA`, and `NA != '+'` is `NA` | filter with `!(ev$Reverse %in% '+')`, which is `FALSE` for `NA` |
