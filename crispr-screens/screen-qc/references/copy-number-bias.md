# Copy-Number Amplicon Bias Diagnostic

## Copy-Number Amplicon Bias Diagnostic

**Goal:** Detect the Aguirre 2016 / Munoz 2016 copy-number artifact where sgRNAs targeting amplified loci appear "essential" purely from DNA-damage burden.

**Approach:** Bin genes by copy number (if known from matched WGS/SNP-array) and check whether mean LFC correlates with CN. Alternatively, count off-target cut sites per sgRNA and check correlation with depletion -- amplified loci share many identical cut sites.

```python
def cn_bias_diagnostic(gene_lfc_df, cn_df):
    '''cn_df: per-gene copy number (from WGS/SNP-array/matched ASCAT).
    Tests whether amplified genes show systematically lower LFC.'''
    merged = gene_lfc_df.merge(cn_df, on='gene')
    bins = pd.qcut(merged['copy_number'], q=5, duplicates='drop')
    bin_lfc = merged.groupby(bins, observed=True)['lfc'].agg(['mean', 'median', 'std', 'count'])
    from scipy.stats import spearmanr, mannwhitneyu
    rho, p = spearmanr(merged['copy_number'], merged['lfc'])
    amplified = merged[merged['copy_number'] > 4]['lfc']
    diploid = merged[merged['copy_number'].between(1.5, 2.5)]['lfc']
    gap, p_gap = np.nan, np.nan
    if len(amplified) >= 3 and len(diploid) >= 3:
        gap = amplified.mean() - diploid.mean()
        p_gap = mannwhitneyu(amplified, diploid, alternative='less').pvalue
    return {'cn_vs_lfc_rho': rho, 'cn_vs_lfc_p': p,
            'n_amplified_genes': len(amplified),
            'amplified_mean_lfc': amplified.mean(),
            'diploid_mean_lfc': diploid.mean(),
            'amplified_vs_diploid_gap': gap,        # negative = amplified genes more depleted
            'p_amplified_more_depleted': p_gap,
            'cn_bias_present': bool((rho < -0.1 and p < 0.01) or (gap < -0.5 and p_gap < 0.01)),
            'per_bin': bin_lfc}
```

**Interpretation: two rules, not one.**

1. **Genome-wide.** Spearman ρ < -0.1 (p < 0.01) between copy number and LFC indicates a broad
   copy-number artifact.
2. **Focal.** Compare `amplified_mean_lfc` with `diploid_mean_lfc` directly. A single amplicon
   covers tens of genes out of ~18,000, so it barely moves ρ: on a realistic 40-gene amplicon the
   genome-wide ρ was only -0.066 while amplified genes averaged LFC -0.877 against -0.019 for
   diploid ones (p = 7e-19). Treat a gap below -0.5 with a significant one-sided test as bias even
   when ρ passes.

Either rule firing means correct before hit calling. When a specific amplicon is suspected, run the
diagnostic again on that region's genes plus a diploid background. Remediation: CRISPRcleanR, CERES
or Chronos (see [[copy-number-correction]], whose `detect_cn_bias()` applies the same two rules)
before hit calling.
