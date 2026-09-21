# CAUSE for CHP-Aware Estimation

## CAUSE for CHP-Aware Estimation

**Goal:** Distinguish causal from shared-factor (correlated horizontal pleiotropy) explanations of an exposure-outcome association.

**Approach:** Fit nuisance parameters (LD pruning + rho_GWAS sample-overlap correction) on a random SNP set; fit the sharing and causal posterior; compare ELPD (expected log predictive density) via Pareto-k smoothed importance sampling.

```r
library(cause)
params <- est_cause_params(dat_cause, variants = pruned_subset_snps)
res_cause <- cause(X = dat_cause, variants = pruned_snps, param_ests = params)
elpd <- summary(res_cause)$tab
```

Full posterior extraction + reporting: examples/cause_analysis.R.

**Interpreting CAUSE output:**

- `q`: posterior CHP fraction; 0 = no CHP, 1 = all instruments operate via the shared factor
- `eta`: shared-factor effect on Y (the "confounder pathway" magnitude)
- `gamma`: posterior causal effect of E on Y after partialling out CHP; report median + 95% credible interval
- `delta_ELPD` (sharing - causal): negative -> causal model preferred; z = delta_elpd / se(delta_elpd); z > 1.96 standard, z > 3.0 stringent; one-sided p reported alongside posterior gamma
- Pareto-k > 0.7 indicates unstable posterior on those points; if more than 10% of points are unstable, treat the posterior as unreliable; remediation: add more SNPs (loosen p-threshold one notch then re-prune in LD) or re-fit excluding flagged outliers

CAUSE requires sumstats from both exposure and outcome GWAS in matched effect-allele coding. The pruning step typically retains 100-5000 signature SNPs at LD r^2 < 0.01 in a 1 Mb window; nuisance estimation should use a larger random SNP subset (`>= 100,000` genome-wide SNPs) to fit rho (sample overlap) stably.
