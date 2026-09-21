# Phenome-Wide Drug-Target MR

Moved verbatim from SKILL.md. Read when scanning one target across many outcomes for on-target adverse effects.

## Phenome-Wide Drug-Target MR

**Goal:** For a single drug target (single protein), test causal effect across hundreds of outcomes to discover on-target adverse effects.

**Approach:** Hold cis-pQTL instrument set fixed; loop outcome over OpenGWAS catalogue or FinnGen DF12; multi-test correct over outcomes.

```r
library(TwoSampleMR); library(ieugwasr)

target_pqtl <- read.table('pcsk9_cis_pqtls.tsv', header = TRUE)
exposure_dat <- format_data(target_pqtl, type = 'exposure',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2', eaf_col = 'EAF', pval_col = 'P')

curated_endpoints <- read.table('finngen_DF12_endpoints.tsv', header = TRUE)  # ~3000 curated endpoints from finngen.fi
outcomes <- available_outcomes()
outcomes_filt <- subset(outcomes, id %in% curated_endpoints$id & sample_size >= 50000 & population == 'European')

results <- lapply(outcomes_filt$id, function(out_id) {
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = out_id)
    if (nrow(outcome_dat) < 2) return(NULL)
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    mr(dat, method_list = 'mr_ivw')
})

results_df <- do.call(rbind, Filter(Negate(is.null), results))
n_tests <- nrow(curated_endpoints)
results_df$p_bonf <- pmin(results_df$pval * n_tests, 1)
top_hits <- subset(results_df, pval < 0.05 / n_tests)   # 0.05 / 3000 = 1.7e-5
```

Document the curated endpoint list (FinnGen DF12, Open Targets curated trait map, or a manuscript-specific phecode hierarchy) in methods. The `outcomes_filt$id[1:200]` pattern is a debug shortcut, not a defensible pheWAS protocol. The PCSK9 -> T2D signal (Schmidt 2017 Lancet Diabetes Endocrinol 5:97) was discovered exactly via curated-endpoint pheWAS: cis-MR of LDL-lowering instruments revealed on-target T2D risk before clinical trials confirmed it. Drug-target pheWAS is the canonical use case.

