## Specificity Test (CellPhoneDB v5)

**Goal:** Get permutation specificity p-values with rigorous multi-subunit complex handling (human).

**Approach:** Run the statistical method on log-normalized counts plus a cell-type meta table; the permutation null shuffles cluster labels, and complexes require all subunits via the limiting (minimum) subunit.

```python
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# threshold=0.1: a gene must be expressed in >=10% of a cluster's cells to count
# iterations=1000: label-permutation null; pvalue=0.05 reports per-pair significance
# debug_seed fixes the permutation RNG for reproducible p-values (default -1 is unseeded),
# but only at threads=1: with threads>1 the pool workers do not replay the RNG stream
# (2 identical seeded runs at threads=4 differed on 82/120,375 significance flags)
# score_interactions=True uses multiprocessing.Pool internally, so the __main__ guard below
# is required on Windows -- without it the call crashes with RuntimeError
def main():
    results = cpdb_statistical_analysis_method.call(
        cpdb_file_path='cellphonedb.zip',          # cellphonedb-data v5 release
        meta_file_path='meta.tsv',                  # barcode -> cell_type
        counts_file_path='counts_normalized.h5ad',  # normalized, NOT scaled
        counts_data='hgnc_symbol',
        threshold=0.1, iterations=1000, pvalue=0.05, debug_seed=1337,
        score_interactions=True, threads=1, output_path='cpdb_out')
    return results
    # DEG-driven escape from one-vs-rest: cpdb_degs_analysis_method.call(..., degs_file_path=...)

if __name__ == '__main__':
    results = main()
```
