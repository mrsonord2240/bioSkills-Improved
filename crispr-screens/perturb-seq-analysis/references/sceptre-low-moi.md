## SCEPTRE for Low-MOI Differential Expression

**Why this matters:** Standard differential-expression tools (DESeq2, MAST) assume Gaussian-mixture distribution and fail at single-cell scale with sparse, zero-inflated data. SCEPTRE (Katsevich Lab, 2021; low-MOI variant Barry 2024 Genome Biol) uses a negative-binomial GLM with conditional resampling:

1. Per gene, fit NB GLM: `log(expr_g) ~ pert_indicator + technical_factors`
2. Compute z-score for the perturbation coefficient
3. Resample the pert_indicator (conditional on counts) 500-1000 times; compute permutation null
4. Get FDR via permutation; not parametric

```r
library(sceptre)

# Input: sce object or sparse matrix + metadata
# Required: gene_expression_matrix, perturbation_indicator (binary per cell per pert),
#           technical_factors (batch, n_genes, etc.)

# For each gene + perturbation pair:
# Current sceptre API is a pipeline of composable steps:
sceptre_object <- import_data(response_matrix, grna_matrix, grna_target_data_frame,
                              moi = 'low', extra_covariates = covariates_df)
sceptre_object <- set_analysis_parameters(sceptre_object, discovery_pairs = pairs_df)
sceptre_object <- assign_grnas(sceptre_object)
sceptre_object <- run_qc(sceptre_object)
sceptre_object <- run_calibration_check(sceptre_object)
sceptre_object <- run_discovery_analysis(sceptre_object)
results <- get_result(sceptre_object, analysis = 'run_discovery_analysis')
# Output: per-gene-per-pert p-value, log-fold-change, FDR
```

**Advantage over MAST:** SCEPTRE's permutation NB GLM is the only method that maintains calibrated FDR in pooled-screen scRNA-seq (Barry 2024 benchmark). MAST and Wilcoxon are over-confident due to data sparsity.
