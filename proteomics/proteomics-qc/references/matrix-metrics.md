# Matrix-level QC metrics: replicate correlation, CV, missingness

Read when computing replicate correlation, a sample-swap check, CV, or a completeness filter, or when choosing how to impute.

## Replicate Correlation on log2

**Goal:** Quantify reproducibility without letting a few abundant proteins fake agreement.

**Approach:** Correlate on log2 intensities (variance-stabilized, high-abundance tail compressed), report within-group pairs, and flag a sample correlating better with another group as a possible swap.

```python
from itertools import combinations

def replicate_correlation(log2_intensities, sample_groups):
    # DESIGN REQUIREMENT: >=2 samples in at least one group. Without it this silently returned
    # an empty table, which reads like "no reproducibility problem" when it means "not measured".
    sizes = sample_groups.value_counts()
    if (sizes >= 2).sum() == 0:
        raise ValueError(f'replicate_correlation needs >=2 samples in a group; sizes are '
                         f'{sizes.to_dict()}. With one run per condition there is no replicate '
                         'reproducibility to measure -- do not report "no outliers found".')
    corr = log2_intensities.corr(method='pearson')  # log2 first: Pearson on raw is a high-abundance artifact
    rows = []
    for group in sample_groups.unique():
        members = sample_groups[sample_groups == group].index
        if len(members) < 2:
            # a row per unchecked group, so the caller sees "not measured" in the data, not only on stdout
            print(f'WARNING: group {group!r} has no within-group pair; its sample is UNCHECKED here.')
            rows.append({'group': group, 's1': members[0], 's2': None, 'r': np.nan, 'status': 'not_measurable_n1'})
        for s1, s2 in combinations(members, 2):
            rows.append({'group': group, 's1': s1, 's2': s2, 'r': corr.loc[s1, s2], 'status': 'measured'})
    return pd.DataFrame(rows)  # filter status == 'measured' before summarizing r

def cross_group_correlation(log2_intensities, sample_groups, top_n=300):
    # Sample-swap check: does a sample match ANOTHER group better than its own? Within-group r alone
    # cannot show it (a swapped pair still correlates 0.93-0.96 with everything). Centre each protein
    # on its mean over samples, on the most variable proteins, so condition -- not shared abundance -- drives r.
    complete = log2_intensities.dropna(how='any')
    top = complete.loc[complete.var(axis=1).nlargest(top_n).index]
    corr = top.sub(top.mean(axis=1), axis=0).corr()
    rows = []
    for s in corr.columns:
        mean_r = {g: corr.loc[s, [x for x in corr.columns if sample_groups[x] == g and x != s]].mean()
                  for g in sample_groups.unique()}
        own = sample_groups[s]
        others = {g: r for g, r in mean_r.items() if pd.notna(r)}
        best = max(others, key=others.get)
        measurable = pd.notna(mean_r[own])
        rows.append({'sample': s, 'own_group': own, 'r_own': mean_r[own], 'best_group': best,
                     'r_best': others[best], 'status': 'measured' if measurable else 'not_measurable_n1',
                     'possible_swap': bool(best != own) if measurable else None})
    return pd.DataFrame(rows)
```

Summarize only `status == 'measured'` rows: a singleton group appears as `not_measurable_n1` with r = NaN, which means UNCHECKED, not clean. A swapped pair still correlates 0.93-0.96 with its own group's other members, so within-group r cannot show a swap; `cross_group_correlation` reports each sample's mean centred r to its own group and to every other group and sets `possible_swap` when another group matches best. A flagged sample is a swap or relabel candidate, not proof: confirm against the sample sheet.

Technical replicates r > 0.98 (instrument noise only); biological r ~ 0.90-0.98 (genuine variance, lower is expected and correct); soft floor r > 0.8 to retain a biological replicate. A Spearman check is a robustness aid only -- ranks discard the magnitude that quant QC cares about.

## Coefficient of Variation on the Linear Scale

**Goal:** Summarize per-condition precision with a number that means what it says.

**Approach:** Compute CV = SD/mean on LINEAR (non-log) intensities; if only logged values exist use the geometric-CV formula. Report the median CV per condition (the per-protein distribution is right-skewed).

```python
def median_cv_linear(linear_intensities, sample_groups):
    # DESIGN REQUIREMENT: >=2 samples per group. A one-sample group yields NaN, and a table of
    # NaNs is not "excellent precision" -- fail loudly instead of returning it.
    sizes = sample_groups.value_counts()
    if (sizes >= 2).sum() == 0:
        raise ValueError(f'median_cv_linear needs >=2 samples in a group; sizes are {sizes.to_dict()}. '
                         'CV is undefined with no replicates -- report "not measurable", not NaN.')
    rows = []
    for group in sample_groups.unique():
        members = sample_groups[sample_groups == group].index
        block = linear_intensities[members]
        per_protein_cv = block.std(axis=1) / block.mean(axis=1)  # base CV formula REQUIRES linear scale
        n1 = len(members) < 2
        if n1:
            print(f'WARNING: group {group!r} has {len(members)} sample(s); its CV is NaN (undefined), not low.')
        rows.append({'group': group, 'median_cv_pct': 100 * per_protein_cv.median(),
                     'status': 'not_measurable_n1' if n1 else 'measured'})
    return pd.DataFrame(rows)

def geometric_cv_from_log(log_intensities):
    sigma = log_intensities.std(axis=1) * np.log(2)  # convert log2 SD to natural-log SD
    return 100 * np.sqrt(np.expm1(sigma ** 2))  # gCV = sqrt(exp(sigma^2) - 1)
```

Applying the base formula to log-transformed data is meaningless (Brenes 2024; see "CV computed on log-transformed data"). A group with one sample returns `status = 'not_measurable_n1'` and NaN: not measurable, not low. State normalization state, transform, and software params or the CV is uninterpretable: DIA-NN "High precision" mode silently median-normalizes, halving median CV vs "High accuracy". Technical median CV < ~10-20%, biological ~20-40%; a LOWER CV is not automatically better (loose FDR or faulty MS1 extraction produce artificially low CVs).

## Missingness Mechanism and Completeness

**Goal:** Decide how to impute by first deciding why values are missing.

**Approach:** Diagnose the missingness profile -- left-tail concentration means MNAR (left-censored, abundance-dependent), all-abundance scatter means MCAR -- and filter on completeness before imputing only the shallow remainder.

```python
def missingness_profile(log2_intensities, n_bins=10):
    present_per_protein = log2_intensities.notna().mean(axis=1)
    mean_abundance = log2_intensities.mean(axis=1)
    abundance_bin = pd.qcut(mean_abundance, n_bins, duplicates='drop')
    # present fraction per mean-abundance bin: rising-with-abundance = MNAR, flat = MCAR
    return present_per_protein.groupby(abundance_bin, observed=True).mean()

def completeness_filter(log2_intensities, sample_groups, min_valid_frac=0.7):
    keep = pd.Series(False, index=log2_intensities.index)
    for group in sample_groups.unique():
        block = log2_intensities[sample_groups[sample_groups == group].index]
        keep |= block.notna().mean(axis=1) >= min_valid_frac  # valid in >=70% of >=1 condition
    return log2_intensities[keep]
```

kNN-imputing a genuinely-absent (MNAR) value invents mid-range abundance and KILLS a real present/absent difference; a left-shifted draw (Perseus down-shifted normal, downshift=1.8 SD below the observed mean, width=0.3 of observed SD) on an MCAR gap FABRICATES a false low and inflates a difference. Match the imputer to the mechanism. The imputation mechanics themselves are quantification.
