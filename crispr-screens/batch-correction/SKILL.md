---
name: bio-crispr-screens-batch-correction
description: Batch effect correction for CRISPR screens covering ComBat empirical-Bayes, RUV, SVA, control-sgRNA normalization, and the model-based alternative of including batch as a covariate in MAGeCK MLE or Chronos. Covers screen-specific batch sources (passage cohort, library lot, infection day, sequencing run, Cas9 lot, FBS lot), PCA + variance-decomposition diagnostic to decide if correction is needed, when correction harms biology by over-correcting condition into batch, limma removeBatchEffect for visualization-only correction, and relationship to multi-condition design matrices. Use when combining screens for joint analysis, when passage cohort confounds biology, when DepMap-style panels need Chronos with batch covariates, when picking ComBat vs RUV, or when correction harms biology and should be replaced with explicit covariate modeling.
tool_type: mixed
primary_tool: pyComBat
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pyComBat 0.3.3+ (epigenelabs/pyComBat), MAGeCK 0.5.9+, R/limma 3.58+, sva 3.50+, RUVSeq 1.36+, pandas 2.2+, numpy 1.26+, scikit-learn 1.4+, scipy 1.12+.

Install: `pip install combat` (provides `combat.pycombat`; the PyPI package named `pycombat` is a different project); `mageck` from bioconda (`conda install -c bioconda mageck`, not on PyPI); R: `BiocManager::install(c('sva', 'RUVSeq', 'limma'))`. Inputs: a count matrix (rows = sgRNA, columns = samples), a metadata table with `batch`, `condition` and `replicate`, and for NTC-anchored normalization a list of non-targeting sgRNAs. Code checked 2026-09-21 on pyComBat (`combat`) 0.3.3, MAGeCK 0.5.9.5, sva 3.54.0, RUVSeq 1.40.0.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show combat`; `from combat.pycombat import pycombat`
- R: `packageVersion('sva')`; `?ComBat`; `packageVersion('RUVSeq')`; `?RUVg`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Batch Correction for CRISPR Screens

**"Correct batch effects in my CRISPR screens"** -> Diagnose the batch source, decide whether to remove via empirical-Bayes (ComBat), explicit covariate modeling (MAGeCK MLE / Chronos design matrix), control-guide-anchored normalization, or unwanted-variation decomposition (RUV, SVA), then apply only the correction that preserves biological condition signal.

- Python: `pyComBat.pycombat` for empirical-Bayes correction
- Python: explicit batch covariates in `mageck mle --design-matrix`
- R: `sva::ComBat`, `RUVSeq::RUVg`, `limma::removeBatchEffect`
- Python: Chronos (`crispr_chronos`) natively handles screen-batch covariates

## Batch Sources in CRISPR Screens

| Source | Mechanism | Detectable by |
|--------|-----------|---------------|
| Library lot | Different aliquots or PCR amplifications | Gini shift; plasmid-pool sequencing |
| Cell passage cohort | Cells passaged through different periods | PCA Day-0 samples clustering by passage |
| Infection day | Lentivirus titer drifts; FBS lot changes | PCA Day-0 samples cluster by day |
| Cas9 enzyme lot | Cas9 expression heterogeneity | PR-AUC drift across screens |
| Sequencing run | Lane bias, flowcell variant, machine | Per-sample read-count distribution |
| FBS / culture lot | Fetal bovine serum lot variations confound proliferation | Day-0 vs endpoint differential not present in vehicle |
| Tissue-prep batch | In-vivo: animal cohort, surgical day, organ-prep tech | In-vivo screens (see [[in-vivo-screens]]) |

**Critical:** Batch effects in CRISPR screens often correlate with biology (e.g., the drug arm was processed in batch 2 because that's when the drug arrived). This confounds correction. Always check for confounding before applying ComBat.

## Batch Effect Decision Tree

| Diagnostic finding | Recommended correction |
|--------------------|------------------------|
| PCA shows samples cluster by condition, not batch | No correction needed; biology dominates |
| PCA PC1 separates batches, PC2 separates conditions | Apply ComBat with condition passed as `mod` |
| Batch fully confounded with condition (e.g. all drug in batch 2, all vehicle in batch 1) | Correction will destroy biology; instead redesign next screen with cross-batch balance OR re-analyze with batch in MAGeCK MLE design matrix |
| Day-0 (pre-perturbation) samples cluster by batch | Strong batch effect; ComBat needed |
| Endpoint samples cluster by batch but not Day-0 | Selection-driven artifact (FBS lot etc); correct or include batch as covariate |
| Replicates within a batch are tight; across-batch much wider | Classic batch effect; ComBat |
| Each replicate scatters randomly across PCs | Sample-level noise; no batch correction will help |
| Cancer-line panel with multiple batches | Use Chronos (built-in batch and CN modeling). Copy-number bias is a separate, batch-like effect per line: apply CN correction (CRISPRcleanR / Chronos) before batch correction |
| Several screens sharing one library | JACKS (joint efficacy across screens) or Chronos |

## Diagnose: PCA + Variance Decomposition

**Goal:** Quantify what fraction of variance is batch vs condition before correcting.

**Approach:** Run PCA on log10(counts+1); fit ANOVA decomposing variance into batch and condition components; report variance explained.

```python
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy import stats

def batch_diagnostic(counts_df, metadata_df, batch_col='batch', condition_col='condition'):
    '''Variance decomposition: report fraction of PC1/PC2 variance attributable to batch vs condition.'''
    log_counts = np.log10(counts_df + 1).T  # samples as rows
    pca = PCA(n_components=5)
    pcs = pca.fit_transform(log_counts)
    out = pd.DataFrame({
        'PC': range(1, 6),
        'var_explained': pca.explained_variance_ratio_,
    })
    pc_df = pd.DataFrame(pcs, columns=[f'PC{i+1}' for i in range(5)], index=counts_df.columns).join(metadata_df)
    for i in range(5):
        pc = pc_df[f'PC{i+1}']
        f_b, p_b = stats.f_oneway(*[pc[pc_df[batch_col] == b] for b in pc_df[batch_col].unique()])
        f_c, p_c = stats.f_oneway(*[pc[pc_df[condition_col] == c] for c in pc_df[condition_col].unique()])
        out.loc[i, 'batch_F'] = f_b
        out.loc[i, 'batch_p'] = p_b
        out.loc[i, 'cond_F'] = f_c
        out.loc[i, 'cond_p'] = p_c
    return out
```

**Interpretation:** If PC1 has batch F-stat > condition F-stat by 10x, batch is dominating and correction is warranted. If condition dominates PC1, no correction needed.

## Related but out of scope

ComBat, RUV and SVA are general methods and the code here would run on bulk RNA-seq or proteomics
matrices, but everything that makes this Skill a *screen* Skill is CRISPR-specific: the NTC-count
rules, the CEGv2 PR-AUC and essential-dropout validation, and the MAGeCK MLE / Chronos integration.
For batch correction outside CRISPR screens, keep the method and replace those checks with the
assay's own.

## ComBat Empirical-Bayes Correction

**Goal:** Remove batch-specific location and scale shifts while preserving biological condition signal.

**Approach:** Log-transform counts, fit ComBat with the condition passed as `mod` (so the model knows which signal to preserve), back-transform. Features ComBat cannot fit are left as raw counts and returned by name.

```python
import numpy as np
import pandas as pd
from combat.pycombat import pycombat

def combat_correct(counts_df, batch_vector, condition_vector=None, verbose=True):
    '''ComBat on log2 counts with optional condition covariate (`mod`).
    Preserves condition signal while removing batch shifts.
    Returns (corrected_counts, uncorrected): the matrix (same rows and columns as the input) and
    the index of features that were left as raw counts.

    pycombat takes the matrix as a DataFrame (rows = features, columns = samples) and both
    `batch` and `mod` as plain lists of labels -- it one-hot-encodes `mod` itself, so passing a
    pre-encoded array fails inside pycombat.

    ComBat divides each feature by its pooled residual variance after batch (and `mod`) are
    regressed out. A feature with none left -- constant in every batch, or fully explained by
    batch + condition (e.g. counts 0,0,1,1 in each batch with condition 0,0,1,1) -- corrupts the
    shared empirical-Bayes prior, with no exception and exit code 0. Depending on the data the
    result is an all-NaN matrix (568,720 of 568,720 values on real TKOv3 counts) or the feature
    silently collapsing to a constant with no NaN at all. The pre-fit filter below is what
    prevents both; the NaN check after the fit is only a backstop. A feature that is constant in
    ONE batch but varies in another is safe and is corrected.
    '''
    data = pd.DataFrame(np.log2(counts_df.values + 1),
                        index=counts_df.index, columns=counts_df.columns)

    # Pre-fit filter: residual sum of squares after regressing out the design pycombat will use.
    design = [pd.get_dummies(pd.Series(list(batch_vector)), dtype=float)]
    if condition_vector is not None:
        design.append(pd.get_dummies(pd.Series(list(condition_vector)), drop_first=True, dtype=float))
    X = pd.concat(design, axis=1).to_numpy()
    resid = data.values - data.values @ (X @ np.linalg.pinv(X))
    usable = pd.Series((resid ** 2).sum(axis=1) > 1e-8, index=data.index)   # not == 0: round-off
    uncorrected = data.index[~usable]
    if verbose and len(uncorrected):
        print(f'ComBat: {len(uncorrected)} features have no variance left after batch and condition; '
              f'returned as raw counts (see the `uncorrected` return value)')

    if condition_vector is not None:
        corrected = pycombat(data[usable], list(batch_vector), mod=list(condition_vector))
    else:
        corrected = pycombat(data[usable], list(batch_vector))

    # Backstop only -- fail loudly rather than pass a dead matrix on.
    if corrected.isna().any().any():
        raise ValueError('ComBat returned NaN values despite the pre-fit filter; do not use the '
                         'output. Inspect the design (batch/condition labels) and the count matrix.')

    out = pd.DataFrame(np.power(2, corrected.values) - 1,
                       index=corrected.index, columns=corrected.columns).clip(lower=0)
    return out.reindex(counts_df.index).fillna(counts_df), uncorrected   # dropped features keep raw counts
```

Run `combat_correct()` only when the diagnostic says batch dominates; on a batch-free design ComBat has nothing to remove and adds noise. Keep the `uncorrected` index and flag those features (or exclude them) in hit calling: they still carry their batch effect. Model them with batch as a covariate instead ("Batch as Explicit Covariate" below).

**Critical caveat:** ComBat assumes batch effects are linear shifts of mean and variance in log space. Non-linear effects (e.g., gene-specific batch sensitivity) remain. Always re-check PCA after correction to confirm batches now overlap.

## RUV (Remove Unwanted Variation)

**Goal:** Identify hidden batch sources via control sgRNAs whose true signal is known.

**Approach:** Designate non-targeting controls as "negative controls" (assumed unchanged); RUV decomposes their variance into unwanted factors, then subtracts these from all data.

```r
library(RUVSeq)
# counts_df: rows = sgRNAs, columns = samples
# RUVg's SeqExpressionSet method takes cIdx as control ROWNAMES (character), not positions;
# which() returns integers and fails S4 dispatch here (the matrix method would accept them).
ntc_rownames <- rownames(counts_df)[rownames(counts_df) %in% ntc_sgrna_names]
stopifnot(length(ntc_rownames) > 0)
seqset <- newSeqExpressionSet(counts = as.matrix(counts_df))
ruv_corrected <- RUVg(seqset, cIdx = ntc_rownames, k = 2)  # k = 2 unwanted factors
# Two outputs. The W factors are what goes into a downstream model (MAGeCK MLE design matrix,
# edgeR/DESeq2 design); normCounts() is the adjusted matrix for PCA and visual checks.
W <- pData(ruv_corrected)          # W_1, W_2: the estimated unwanted factors, one column per k
corrected_counts <- normCounts(ruv_corrected)
```

**When to use:** RUV preferred over ComBat when batches are not annotated (e.g., unknown technical confounders). Worse than ComBat when batch is known and well-annotated; ComBat is more direct.

## SVA (Surrogate Variable Analysis)

**Goal:** Estimate unknown latent factors that may confound the screen.

**Approach:** SVA computes surrogate variables that capture variance not explained by known biological factors; these can then be added to the MAGeCK MLE design matrix as covariates.

```r
library(sva)
# counts_df: rows = sgRNAs, columns = samples
mod <- model.matrix(~ condition, data = metadata)
mod0 <- model.matrix(~ 1, data = metadata)
sv_obj <- sva(as.matrix(counts_df), mod, mod0)
n_sv <- sv_obj$n.sv  # number of surrogate variables
# Add to design matrix for MAGeCK MLE
design_mat <- cbind(mod, sv_obj$sv)
```

**Use case:** When the screen has clear biological signal (e.g. essentiality recovery passes) but small effect sizes are hidden by noise; SVA-discovered latent factors as covariates can recover them.

## Batch as Explicit Covariate (Preferred for MAGeCK MLE / Chronos)

**Goal:** Model batch and biology in the same regression instead of pre-correcting.

**Approach:** Add batch indicator columns to the MLE design matrix. The fitted beta for condition is the effect after accounting for batch; no pre-correction needed.

```bash
# Design matrix for a screen with 2 batches and 2 conditions
cat > design.txt <<EOF
Samples         baseline    batch2    treatment
Veh_b1_r1       1           0         0
Veh_b1_r2       1           0         0
Drug_b1_r1      1           0         1
Drug_b1_r2      1           0         1
Veh_b2_r1       1           1         0
Veh_b2_r2       1           1         0
Drug_b2_r1      1           1         1
Drug_b2_r2      1           1         1
EOF

mageck mle \
    --count-table counts.txt \
    --design-matrix design.txt \
    --permutation-round 10 \
    --output-prefix batch_aware_mle
```

**Runtime:** at genome scale (~18,000 genes) `mageck mle`'s variance-model permutation is a multi-hour
job by design (a long run is not a hang). For a fast sanity check, run a gene subset or lower
`--permutation-round`; use the full run for the reported result.

**Reproducibility:** the beta estimates are deterministic, but `mageck mle`'s significance comes from
a permutation procedure that is not seeded, so p-values and FDRs move slightly between reruns on
identical input. Fix `--permutation-round` (higher = more stable, linearly slower) and report the
value, or treat borderline FDRs as borderline.

**Why this is preferred:** ComBat shifts counts before testing; the MLE-with-covariates approach correctly propagates uncertainty from the batch term into the condition beta's standard error. ComBat-then-test pretends the corrected counts are noise-free, biasing FDR.

## Control-Sgrna Anchored Normalization

**Goal:** Use non-targeting controls as the per-sample reference so batch shifts cancel.

**Approach:** Scale each sample so its NTC sgRNAs have a constant median. Subsequent fold changes are relative to NTCs in each sample, automatically batch-controlling.

```python
def ntc_anchored_normalize(counts_df, ntc_sgrna_names, target_median=1000):
    '''Scale each sample so its NTC median is target_median. Subsequent LFC is NTC-anchored.'''
    is_ntc = counts_df.index.isin(ntc_sgrna_names)
    ntc_medians = counts_df.loc[is_ntc].median(axis=0)
    scale_factors = target_median / ntc_medians.replace(0, np.nan)
    return counts_df * scale_factors, scale_factors
```

**Critical:** Requires ≥500 NTCs in the library (see [[library-design]]). With fewer, the NTC median is unstable and amplifies noise rather than removing batch; fall back to median normalization.

## When NOT to Correct

| Situation | Why correction hurts |
|-----------|----------------------|
| Batch is fully confounded with condition | Correction destroys biology along with batch; redesign or accept |
| Batch effect is smaller than between-replicate noise | Correction adds noise without removing meaningful variance |
| Replicates already correlate >0.95 within and across batches | No batch effect to correct |
| Single-screen analysis | No "batch" to correct; only replicate noise |
| Per-batch sample size <3 | Cannot estimate batch shift reliably; correction is harmful |

## Failure Modes

### ComBat eliminates biological signal

**Trigger:** Batch is correlated with condition (e.g., all drug-arm samples were processed week 2; all vehicle-arm samples week 1).
**Mechanism:** ComBat without a `mod` covariate treats condition variance as batch variance; corrects it away.
**Symptom:** PR-AUC against CEGv2 drops after ComBat correction.
**Fix:** Always supply `mod` (the condition labels); verify by comparing PR-AUC before and after.

### RUV adds noise instead of removing it

**Trigger:** k (number of unwanted factors) set too high.
**Mechanism:** RUV's least-squares decomposition over-fits; "removed" variance includes biology.
**Symptom:** Hits decrease and replicate Pearson drops after correction.
**Fix:** Choose k via cross-validation; default k=1 or 2 for most screens.

### Batch-aware MLE collinear design matrix

**Trigger:** Adding a batch indicator that is fully collinear with another design column (e.g., all of batch 2 is also Day 21).
**Mechanism:** MLE design matrix is singular; betas not estimable.
**Symptom:** MAGeCK MLE errors out or produces NaN betas.
**Fix:** Drop the collinear column; re-design experiment with cross-batch balance.

### ComBat after RUV double-corrects

**Trigger:** Applying multiple corrections sequentially.
**Mechanism:** Both methods remove variance; sequential application removes biology twice.
**Symptom:** All signal gone; counts look uniformly noisy.
**Fix:** Pick one method based on diagnostic; never combine.

### Per-batch sample size too small

**Trigger:** 2 replicates per batch with 3 batches; ComBat estimates batch shift from 2 samples.
**Mechanism:** Insufficient data to estimate batch parameters; high-variance estimates.
**Symptom:** Correction makes some batches worse than uncorrected.
**Fix:** Need ≥3 (preferably 4-6) samples per batch; below this, use covariate modeling instead.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| PC1 batch F vs condition F | F_batch > 10x F_cond -> apply correction | Standard variance-decomposition diagnostic |
| ComBat min samples per batch | ≥3, ideally 4-6 | Empirical Bayes prior estimation |
| RUV `k` (unwanted factors) | k=1 default; k=2 if multiple known batch sources | Risso 2014; cross-validate |
| NTCs needed for NTC-anchored norm | ≥500 in library | Stable median |
| Post-correction PCA check | Batches must overlap in PC1/PC2 plot | Visual sanity check |
| Post-correction PR-AUC | Should be same or higher than pre | If lower, correction destroyed biology |

## Validation Checklist

After applying correction:

- [ ] Corrected matrix has no NaN/Inf values (ComBat can return all-NaN with exit code 0)
- [ ] Features in `uncorrected` are flagged or excluded in hit calling
- [ ] PCA: batches now overlap (visual)
- [ ] Within-batch Pearson preserved (should be unchanged)
- [ ] Across-batch Pearson improved
- [ ] CEGv2 PR-AUC preserved or higher
- [ ] NTC distribution stable across batches
- [ ] No new outlier samples introduced
- [ ] Hit list compared with the uncorrected hit list; every difference explained by the batch effect, not lost biology

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| PR-AUC drops after ComBat | Batch confounded with condition | Add `mod` covariate; or redesign |
| `combat_correct()` reports N features uncorrected | No variance left after batch and condition (e.g. all-zero guides), so ComBat cannot fit them | Expected; flag or exclude the returned `uncorrected` index in hit calling, or model batch as a covariate |
| MAGeCK MLE NaN beta after adding batch column | Collinear design matrix | Drop collinear column |
| Replicates still cluster by batch after RUV | k too low | Increase k; cross-validate |
| Replicates lose internal cohesion after correction | Over-correction | Reduce k or revert |
| NTC-anchored norm worse than median | Too few NTCs | Use median; add NTCs to next library |
| Sequencing-run-level batch survives ComBat | Non-linear sequencing effect | Pre-normalize with `mageck count --norm-method control` first |

## References

- Johnson WE et al. 2007. *Biostatistics* 8:118. Original ComBat algorithm.
- Leek JT et al. 2012. *Bioinformatics* 28:882. SVA package.
- Risso D et al. 2014. *Nat Biotechnol* 32:896. RUVSeq.
- Pacini C et al. 2021. *Nat Commun* 12:1661. Integrated cross-study dependencies; cross-screen batch-effect correction.
- Vinceti A et al. 2024. *Genome Biol* 25:192. Benchmark of methods for correcting biases in CRISPR-Cas9 screening data.

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK MLE with explicit batch covariates
- crispr-screens/screen-qc - Pre-correction PCA diagnostic
- crispr-screens/copy-number-correction - Chronos handles batch + CN jointly
- crispr-screens/library-design - NTC composition for NTC-anchored normalization
- crispr-screens/jacks-analysis - Joint analysis across batches with shared efficacy
- crispr-screens/hit-calling - Post-correction hit calling
- crispr-screens/in-vivo-screens - In-vivo-specific batch sources (animal cohort, tissue prep)
