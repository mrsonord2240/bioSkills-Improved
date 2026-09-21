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
