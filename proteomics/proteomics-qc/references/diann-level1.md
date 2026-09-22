# Level-1 run metrics from a DIA-NN report

Read when a DIA-NN report is the only instrument-level evidence (RT fit, peak width).

## Level-1 Run Metrics From a DIA-NN Report

**Goal:** Read retention-time fit and peak width per run from the columns a DIA-NN report already carries.

**Approach:** Filter to Global.Q.Value <= 0.01, then per run correlate `RT` with `Predicted.RT` (R^2) and take the median `FWHM` (minutes in DIA-NN 2.x) and `Quantity.Quality`. With no rolling baseline yet, the across-run median stands in for it: a run whose FWHM is >25% above that median (the Thresholds row's 20-30% alarm) or whose RT fit R^2 is below 0.99 is flagged.

```python
def diann_level1(report):
    # report: DIA-NN report.tsv/parquet as a DataFrame (columns Run, RT, Predicted.RT, FWHM, Quantity.Quality, Global.Q.Value)
    rep = report[report['Global.Q.Value'] <= 0.01].dropna(subset=['RT', 'Predicted.RT', 'FWHM'])
    rows = []
    for run, d in rep.groupby('Run'):
        rows.append({'run': run, 'n_precursors': len(d),
                     'rt_fit_r2': np.corrcoef(d['RT'], d['Predicted.RT'])[0, 1] ** 2,
                     'median_fwhm': d['FWHM'].median(),
                     'median_quantity_quality': d['Quantity.Quality'].median()})  # reported only: no accepted cutoff
    out = pd.DataFrame(rows).set_index('run')
    out['fwhm_vs_median'] = out['median_fwhm'] / out['median_fwhm'].median()
    out['flag'] = (out['rt_fit_r2'] < 0.99) | (out['fwhm_vs_median'] > 1.25)
    return out
```

Three runs are a weak baseline: one broad run among three shifts the median. Trend these numbers over the queue (Levey-Jennings) before acting on a single flag.
