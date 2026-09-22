# PCA, batch detection and exclusion rules

Read when checking for batch effects, deciding whether to exclude a sample, or before asking Sam about exclusions (the stop conditions are in the closing paragraph).

## PCA and Batch Detection

**Goal:** See whether the dominant variance is biology or batch, and flag outlier samples.

**Approach:** On the normalized survivors, run PCA, color by condition and by batch, and test whether top PCs associate with batch.

```python
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import f_oneway

def pca_batch_check(normalized_log2, sample_info, batch_col='batch'):
    # sample_info must be indexed by sample name, e.g. pd.read_csv(...).set_index('sample')
    if set(sample_info.index) != set(normalized_log2.columns):
        raise ValueError('sample_info index must equal the matrix column names (set_index on the sample column)')
    # complete cases only: a row-median fill pulls high-missing (failed) samples to the centre of the PCA
    complete = normalized_log2.dropna(how='any')
    n_samples = complete.shape[1]
    if n_samples < 3 or len(complete) < n_samples:
        raise ValueError(f'too few samples ({n_samples}) or complete proteins ({len(complete)}) for PCA')
    n_pc = min(5, n_samples - 1)
    scaled = StandardScaler().fit_transform(complete.T)
    # 'full' is exact and cheap at QC sizes; the default 'auto' is randomized on a wide matrix (see Common Errors)
    pcs = PCA(n_components=n_pc, svd_solver='full', random_state=0).fit(scaled)
    coords = pd.DataFrame(pcs.transform(scaled), columns=[f'PC{i+1}' for i in range(n_pc)],
                          index=complete.columns).join(sample_info)
    print(f'PCA on {len(complete)} complete proteins of {len(normalized_log2)}')
    tests = []
    for pc in coords.columns[:min(3, n_pc)]:
        groups = [coords[coords[batch_col] == b][pc] for b in coords[batch_col].unique()]
        # a level with one sample makes f_oneway raise 'At least two samples are required; got 1'
        if len(groups) < 2 or min(len(g) for g in groups) < 2:
            print(f'{pc} ~ {batch_col}: NOT TESTABLE, level sizes {[len(g) for g in groups]} '
                  f'(need >=2 levels with >=2 samples each) -- this is not evidence of no batch effect')
            tests.append({'pc': pc, 'p': np.nan, 'status': 'not_testable'})
            continue
        _, p = f_oneway(*groups)
        print(f'{pc} ~ {batch_col}: p={p:.4f}')
        tests.append({'pc': pc, 'p': p, 'status': 'tested'})
    return coords, pcs.explained_variance_ratio_, pd.DataFrame(tests)  # tests.status: tested / not_testable
```

The third return value `tests` carries the per-PC status (`tested` / `not_testable`); a `not_testable` row is never evidence of no batch effect. A sample isolated from its group is a removal/re-run candidate, but a high-missing sample is judged by `raw_sample_qc`, not by PCA. If batch is PC1, keep batch in the design matrix for the differential test (preferred when batch and condition are balanced), and use `limma::removeBatchEffect` (or ComBat) only on the matrix used for PCA/plots to re-inspect biology; do not test on a batch-corrected matrix and also model batch. If batch is FULLY confounded with condition (every batch level holds exactly one condition) nothing can be corrected: batch and condition are the same variable, and removing one removes the other -- on a fully confounded synthetic set the mean |log2FC| of 104 truly-changed proteins went from 1.55 to 0.00 after batch removal. Report the design as non-identifiable and stop; do not correct, and do not test. Document and justify every exclusion, and re-run the downstream check with and without borderline samples (`references/qc-report-template.md` has the exclusion log and the with/without table). Stop and ask before excluding samples, when n < 5 per group makes PCA unstable, or when no un-normalized column is available for the loading check. Visualization of the projection routes to data-visualization/dimensionality-reduction-plots.
