# HIMA High-Dimensional EWAS Mediation

### HIMA covariate or data.pheno error

**Trigger:** `data.pheno` contains factor/character columns with NA, or formula references columns missing from `data.pheno`.

**Mechanism:** HIMA v2.3+ uses a formula interface and constructs the design matrix internally from `data.pheno` via `hima_dblasso`'s `process_var()`, which requires every covariate to already be numeric or dummy-coded -- **not** an R `factor`. Passing a `factor` column throws `Non-numeric variable(s) detected... Please convert all factor/character variables to numeric or dummy variables` (verified on HIMA 2.3.4, 2026-09-17); separately, missing values or unparseable formulas surface as cryptic `glmnet`/`storage.mode` errors.

**Symptom:** Pipeline fails inside `hima()` with a non-obvious `storage.mode`/`model.matrix` error, or -- if `factor()` was tried as the fix -- with `process_var()`'s "Non-numeric variable(s) detected" error.

**Fix:** Pre-clean `data.pheno` (drop NA rows for the variables in the formula; ensure all RHS variables exist as columns) and **dummy-code categorical covariates with `model.matrix()`, not `factor()`** -- HIMA 2.3.4 rejects factor columns outright. Example (verified: recovers planted mediators with a 3-level `batch` covariate):
```r
dat <- na.omit(dat[, c('outcome', 'exposure', 'age', 'sex', 'batch', 'pc1', 'pc2')])
batch_dummy <- model.matrix(~ batch, data=dat)[, -1, drop=FALSE]  # drop the intercept column
dat <- cbind(dat[, setdiff(names(dat), 'batch')], batch_dummy)     # e.g. adds batchB, batchC
result <- hima(as.formula(paste('outcome ~ exposure + age + sex +',
                                 paste(colnames(batch_dummy), collapse=' + '), '+ pc1 + pc2')),
               data.pheno=dat, data.M=M_matrix, mediator.type='gaussian')
```

### HIMA mediator-type vs outcome-type mismatch

**Trigger:** Survival outcome (`Surv()` on LHS) with `mediator.type='compositional'` chosen against text mediator panel; or count mediators passed as `'gaussian'`.

**Mechanism:** HIMA v2.3+ auto-detects outcome family from the LHS of `formula` (continuous, binary, survival, count); `mediator.type` is set for the mediator data only (`'gaussian'`, `'negbin'`, `'compositional'`). Mismatching mediator-type to the actual mediator distribution biases the screening step.

**Symptom:** Hazard / rate ratios for indirect effects look implausible; many "significant" mediators fail replication.

**Fix:** Set `mediator.type='gaussian'` for continuous (e.g., methylation beta, log-CPM expression), `'negbin'` for raw count (RNA-seq), `'compositional'` for relative-abundance microbiome. Verify by `?hima` in the installed version since the catalogue of mediator types has expanded across releases.

Cell composition is the canonical unmeasured confounder in EWAS mediation: include estimated cell proportions (Houseman or reference-free RPC method) as covariates -- in the formula interface, add them as RHS terms alongside `age + sex + ...`; in `hima_classic()`, pass them in both `COV.XM` and `COV.MY`.

### High-Dimensional EWAS Mediation (HIMA2)

**Goal:** Among thousands of candidate CpG mediators, identify those mediating an exposure-outcome effect with FDR control.

**Approach:** HIMA v2.3+ uses a formula interface and auto-detects outcome family (Gaussian / binomial / Cox / Poisson). Screening + MCP/DBlasso penalisation + joint significance with BH; the `sigcut` argument controls the FDR threshold (default 0.05).

```r
library(HIMA)

dat <- na.omit(dat[, c('outcome', 'exposure', 'age', 'sex', 'cell_pc1', 'cell_pc2')])
M_matrix <- as.matrix(beta_values)

result <- hima(
  outcome ~ exposure + age + sex + cell_pc1 + cell_pc2,
  data.pheno=dat,
  data.M=M_matrix,
  mediator.type='gaussian',          # 'negbin' for count, 'compositional' for microbiome
  penalty='DBlasso',                  # default; alternatives 'MCP', 'SCAD', 'lasso'
  scale=TRUE,
  sigcut=0.05,
  parallel=TRUE, ncore=8, verbose=TRUE
)
# result is a LIST of class "hima" ($ID, $alpha, $beta, `$alpha*beta`, $rimp, `$p-value`),
# NOT a data.frame -- nrow(result) and rownames(result) both return NULL rather than
# erroring (verified on HIMA 2.3.4). Use result$ID and length(result$ID):
sig_mediators <- result$ID
n_sig <- length(sig_mediators)
```

For survival outcomes wrap the LHS as `Surv(time, status)`; HIMA auto-routes to Cox. The old `hima_classic()` (Zhang 2016 original) is still exported but screens by beta only and misses mediators with strong alpha + weak beta -- prefer the wrapper `hima()` unless reproducing a 2016-2021 paper.

For highly-correlated mediators (CpG-island clusters, gene-module co-expression): HIMA uses joint significance with BH-FDR on max(p_alpha, p_beta) and handles correlation only weakly. Within `hima()`, set `penalty='MCP'` for stronger correlation handling; alternatively pre-reduce the mediator panel by principal components or by clustering correlated mediators and screening the cluster centroid (VanderWeele & Vansteelandt 2014 Epidemiol Methods 2:95).
