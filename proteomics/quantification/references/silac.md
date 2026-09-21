## SILAC Quantification

**Goal:** Compute heavy/light ratios while preserving on/off biology and flagging label artifacts.

**Approach:** A protein present only in the heavy channel may be the interesting biology (or a detection-limit dropout), so do not discard it -- but keep it in a presence flag, not as +/-Inf inside the ratio matrix, where it turns pandas SDs into NaN and stops `eBayes(trend=TRUE, robust=TRUE)`. Verify labeling efficiency (>=95%, target 97-98%) on a heavy-only pilot and assess Arg->Pro conversion before trusting any ratio -- both bias every ratio in the same direction, so neither shows up as extra scatter.

### Check labeling efficiency and Arg->Pro first
```python
import numpy as np, pandas as pd

def silac_labeling_efficiency(pilot, heavy='Intensity H', light='Intensity L'):
    '''Incorporation on a HEAVY-ONLY pilot (cells grown in heavy medium, NOTHING mixed in): every
    light ion there is unlabeled protein. pilot = the pilot's peptide table (MaxQuant evidence.txt).'''
    h = pd.to_numeric(pilot[heavy], errors='coerce').fillna(0.0)
    l = pd.to_numeric(pilot[light], errors='coerce').fillna(0.0)
    ok = (h + l) > 0
    if not ok.any():
        raise ValueError('no peptide has signal in either channel -- wrong columns or wrong file')
    per_pep = h[ok] / (h[ok] + l[ok])
    eff = float(h[ok].sum() / (h[ok] + l[ok]).sum())    # intensity-weighted = the number to report
    # A 1:1 forward mix of these cells does NOT read log2 H/L = 0: the unincorporated (1 - eff) of
    # the heavy sample is counted in the LIGHT channel, so H/L = eff / (2 - eff) -- -0.20 at 93%.
    return {'n_peptides': int(ok.sum()), 'incorporation': round(eff, 4),
            'median_peptide_incorporation': round(float(per_pep.median()), 4),
            'peptides_below_95pct': int((per_pep < 0.95).sum()),
            'expected_log2_HL_bias_at_1to1': round(float(np.log2(eff / (2 - eff))), 4),
            'pass_95pct': bool(eff >= 0.95)}

def arg_to_pro_shift(peptides, seq='Sequence', ratio='Ratio H/L'):
    '''Arg->Pro drains the heavy channel once per proline, so log2 H/L falls with PROLINE COUNT.
    The dose slope is what makes this specific -- a flat offset is incomplete labeling instead.
    Direct route: re-search the pilot with Pro6 variable and take I(Pro6)/(I(Pro6)+I(Pro0)).'''
    r = np.log2(pd.to_numeric(peptides[ratio], errors='coerce'))
    npro = peptides[seq].astype(str).str.count('P')
    ok = np.isfinite(r)
    r, npro = r[ok], npro[ok]
    slope = float(np.polyfit(npro, r, 1)[0]) if npro.nunique() > 1 else float('nan')
    return {'n_pro_free': int((npro == 0).sum()), 'n_pro_containing': int((npro > 0).sum()),
            'median_log2_HL_pro_free': round(float(r[npro == 0].median()), 4),
            'log2_HL_slope_per_proline': round(slope, 4),
            'conversion_per_proline': round(float(1 - 2 ** slope), 4)}
```

Seeded pilot recovering the planted 0.93 incorporation and 0.08 conversion: `examples/lfq_normalization.py`.

### Ratios that keep on/off biology
```python
import numpy as np

# Arg10/Lys8 is the common pairing (avoids overlap with the +6 isotope envelope)
SILAC_SHIFTS = {'Arg10': 10.008269, 'Lys8': 8.014199, 'Arg6': 6.020129, 'Lys6': 6.020129}

def silac_log2_ratio(heavy, light):
    '''Vectorized over arrays/Series: log2 H/L (NaN unless both channels quantified) plus a presence flag.'''
    heavy, light = np.asarray(heavy, dtype=float), np.asarray(light, dtype=float)
    h, l = heavy > 0, light > 0    # NaN compares False
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(h & l, np.log2(heavy / light), np.nan)
    presence = np.select([h & l, h, l], ['both', 'H-only', 'L-only'], default='none')    # report H-only/L-only separately
    return ratio, presence
```
