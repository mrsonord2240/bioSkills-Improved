# q-values from a results table

Read when you hold a PSM table from any engine (concatenated search) or separate target and decoy tables and need q-values. `examples/separate_search_fdr.py` ships the separate-search code as a script.

### FDR from a Results Table (concatenated competition, made explicit)

**Goal:** Compute q-values from any engine's PSM table when the search was a single concatenated target-decoy search.

**Approach:** Keep the best hit per spectrum, rank by score, walk down accumulating target and decoy counts, FDR = (decoys + 1)/targets, then take the running minimum from the bottom to get monotone q-values. `score` must be higher-is-better, and the decoy prefix must match the engine's (Sage and FragPipe write lowercase `rev_`); a table with no recognised decoys must stop, not pass every PSM. This form is correct ONLY for concatenated competition; separate searches need pi0 * decoys/targets (Kall et al. 2008; pi0 = 1 is the conservative default) or the mix-max estimator (Keich, Kertesz-Farkas & Noble 2015; Percolator's default for separate-search input).

```python
import pandas as pd

DECOY_PREFIXES = ('decoy_', 'rev_', 'xxx_')   # compared lower-cased: DECOY_, REV_ / REV__ (MaxQuant), rev_ (Sage, FragPipe), XXX_

psms = pd.read_csv('search_results.tsv', sep='\t')   # map engine columns to 'scan', 'score', 'protein'
# score must be HIGHER-is-better: Sage sage_discriminant_score or Comet xcorr as is;
# E-values (Comet e-value, MS-GF+ SpecEValue) as -log10(E-value)
psms['is_decoy'] = psms['protein'].str.lower().str.startswith(DECOY_PREFIXES)
if not psms['is_decoy'].any():
    raise ValueError(f'no decoy PSMs recognised in {len(psms)} rows: check the decoy prefix, or the table was already decoy-filtered '
                     f'(with no decoys the smallest reachable q is 1/{len(psms)} = {1/len(psms):.3f})')
# one best hit per spectrum (Comet .txt writes 5 rows per scan by default)
psms = psms.sort_values('score', ascending=False).drop_duplicates('scan').reset_index(drop=True)

# concatenated target-decoy competition: each decoy above threshold estimates one false target
targets = (~psms['is_decoy']).cumsum()
decoys = psms['is_decoy'].cumsum()
psms['fdr'] = (decoys + 1) / targets.clip(lower=1)   # +1: zero decoys is not zero FDR (OpenMS conservative default)
psms['qvalue'] = psms['fdr'][::-1].cummin()[::-1]   # running min from the bottom -> monotone q-values

kept = psms[(psms['qvalue'] <= 0.01) & (~psms['is_decoy'])]   # 1% list-level FDR
```

### FDR from SEPARATE Target and Decoy Searches (pi0 * D / T)

**Goal:** Get a valid 1% list out of two result tables produced by searching the same spectra against a target DB and a decoy DB independently.

**Approach:** No competition resolved which hit wins, so the decoy count estimates the number of incorrect TARGETS directly, scaled by pi0, the proportion of target PSMs that are incorrect (Kall et al. 2008). pi0 = 1 is always valid and conservative; the median-decoy estimate (twice the fraction of target scores below the median decoy score) recovers the identifications that pi0 = 1 throws away, at the cost of estimating a nuisance parameter. Do NOT run the concatenated snippet above on the two tables merged. `examples/separate_search_fdr.py` ships this as a script.

```python
import numpy as np

def estimate_pi0(target_scores, decoy_scores):        # Kall et al. 2008
    median_decoy = np.median(decoy_scores)
    return min(1.0, 2.0 * np.mean(np.asarray(target_scores) < median_decoy))

def separate_search_qvalues(targets, decoys, pi0=None):
    # targets, decoys: DataFrames with 'scan' and 'score', ONE row per spectrum each
    if pi0 is None:
        pi0 = estimate_pi0(targets['score'].to_numpy(), decoys['score'].to_numpy())
    t = targets.sort_values('score', ascending=False).reset_index(drop=True)
    decoy_sorted = np.sort(decoys['score'].to_numpy())
    n_decoy_above = len(decoy_sorted) - np.searchsorted(decoy_sorted, t['score'].to_numpy(), side='left')
    t['fdr'] = np.minimum(1.0, pi0 * n_decoy_above / np.arange(1, len(t) + 1))
    t['qvalue'] = t['fdr'][::-1].cummin()[::-1]
    return t, pi0
```

On the synthetic separate-search pair with ground truth (12,000 spectra each): pi0-hat = 0.611 keeps 2,888 PSMs at a true FDP of 1.04%, pi0 = 1 keeps 2,588 at 0.62%, and Elias-Gygi's 2d/(t+d) misapplied here keeps only 2,139 at 0.33%. The alternative is to hand both tables to Percolator and let mix-max do it: that is Percolator's default for separate-search input, and it is a calibrated-score procedure, not the same arithmetic.
