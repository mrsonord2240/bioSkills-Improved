# Second-Best sgRNA Rule and Custom z-score Calling

Read when filtering a hit list for single-guide outliers, or when neither MAGeCK nor BAGEL2 fits the design.

## Second-Best sgRNA Conservative Rule

**Goal:** Reduce false positives from single outlier sgRNAs by requiring the second-most-extreme guide per gene to also be a hit.

**Approach:** For each gene, sort sgRNAs by LFC; require the second-best LFC to exceed a threshold. Rejects genes that depend on one extreme guide.

```python
def second_best_lfc(sgrna_lfc_df, genes_series, direction='neg'):
    '''Return per-gene LFC of the second-best sgRNA in the direction of interest,
    and flag genes with fewer than 2 sgRNAs. For dropout (direction="neg"),
    second-most-negative LFC. A gene with only one sgRNA has no second guide to
    check at all -- return NaN and single_guide=True for it rather than silently
    falling back to the lone guide's own LFC, which would read as "passing" the
    rule with no corroborating guide involved.'''
    results = []
    for gene in genes_series.unique():
        gene_lfc = sgrna_lfc_df[genes_series == gene].sort_values()
        n = len(gene_lfc)
        if n >= 2:
            second = gene_lfc.iloc[1] if direction == 'neg' else gene_lfc.iloc[-2]
            single = False
        else:
            second = float('nan')
            single = True
        results.append({'gene': gene, 'second_best_lfc': second, 'single_guide': single})
    return pd.DataFrame(results)
```

**Rule:** A high-confidence hit has second-best LFC also passing the threshold. A guide-of-one hit has only one extreme guide and should be flagged for orthogonal validation. This rule predates JACKS and is implicit in MAGeCK RRA but explicit elsewhere. Genes with `single_guide=True` (fewer than 2 sgRNAs in the library) have no second guide to check by construction -- always send these to orthogonal validation rather than treating a NaN second-best LFC as a pass.

## Custom z-score Hit Calling (when standard tools don't fit)

**Goal:** Compute gene-level z-scores when neither MAGeCK nor BAGEL2 fits the experimental design.

**Approach:** RPM-normalize, compute per-sgRNA log2 fold-changes, aggregate to gene level, derive z-score from the null distribution of non-targeting controls (cleanest) or all genes (assumes <40% changing), apply BH correction.

```python
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def custom_zscore_hit_calling(counts_df, ctrl_cols, treat_cols, genes_series, ntc_genes=None):
    '''Z-score gene-level hit calling. If ntc_genes provided, null derived from NTCs only;
    otherwise from all genes (assumes <40% changing).'''
    def rpm(df):
        return df.div(df.sum(axis=0), axis=1) * 1e6
    ctrl_rpm = rpm(counts_df[ctrl_cols])
    treat_rpm = rpm(counts_df[treat_cols])
    lfc_per_sgrna = np.log2((treat_rpm.mean(axis=1) + 1) / (ctrl_rpm.mean(axis=1) + 1))
    gene_lfc = pd.DataFrame({'gene': genes_series, 'lfc': lfc_per_sgrna}).groupby('gene')['lfc'].agg(['mean', 'std', 'count'])
    gene_lfc.columns = ['mean_lfc', 'std_lfc', 'n_sgrnas']
    if ntc_genes is not None:
        null = gene_lfc.loc[gene_lfc.index.isin(ntc_genes), 'mean_lfc']
        null_mean, null_std = null.median(), null.std()
    else:
        null_mean = gene_lfc['mean_lfc'].median()
        null_std = gene_lfc['mean_lfc'].std()
    gene_lfc['z'] = (gene_lfc['mean_lfc'] - null_mean) / null_std
    gene_lfc['p'] = 2 * stats.norm.sf(np.abs(gene_lfc['z']))
    gene_lfc['fdr'] = multipletests(gene_lfc['p'], method='fdr_bh')[1]
    return gene_lfc.sort_values('z')
```

