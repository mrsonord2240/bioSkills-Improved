## AP-MS / Affinity-Enrichment Scoring

**Goal:** Rank prey in a pulldown by enrichment over NEGATIVE-CONTROL IPs, not by abundance or by ratio to the input lysate.

**Approach:** A pulldown is deliberately non-representative, so no data-internal normalization (median, sample-loading, IRS) applies, and the input lysate is not a control -- sticky background (ribosome, chaperones, tubulin, keratin) binds the beads in the pulldown and is diluted in the input, so "top N over input" returns background. The control IP is the only thing that separates a bead binder from an interactor. Require reproducible detection across bait replicates and enrichment over control, and floor absent controls at the run's detection limit so bait-only prey score finitely instead of `+Inf`. The default `min_bait_reps` is every bait replicate (reproducibility first), which trades sensitivity for specificity: at 3 replicates a true interactor missing from one by stochastic dropout is dropped (an enrichment of 5.4 was, on the audit fixture). `min_bait_reps=2` is the usual compromise; run both and report how many prey the looser setting adds so the choice is visible.

```python
import numpy as np, pandas as pd

def score_vs_control_ips(ip, bait_cols, ctrl_cols, fc_cutoff=2.0, min_bait_reps=None):
    '''ip: prey x replicate RAW intensities or spectral counts; 0/NaN = not detected. Prey never
    seen in any bait IP get a NaN enrichment (nothing was measured) and are never called.'''
    L = np.log2(ip[list(bait_cols) + list(ctrl_cols)].replace(0, np.nan).astype(float))
    if min_bait_reps is None:
        min_bait_reps = len(bait_cols)          # default: every bait replicate, reproducibility first
    n_bait, n_ctrl = L[bait_cols].notna().sum(axis=1), L[ctrl_cols].notna().sum(axis=1)
    # a control IP that produced no data is skipped by every mean/min below, so say so out loud
    dead = [c for c in ctrl_cols if L[c].notna().sum() == 0]
    if dead:
        print(f'AP-MS: control runs with no data, excluded: {dead}')
    if len(dead) == len(ctrl_cols):
        raise ValueError('every control IP is empty -- nothing to score against')
    Lf = L.copy()
    Lf[ctrl_cols] = Lf[ctrl_cols].fillna(L[ctrl_cols].min())      # per-control-run detection floor
    enrich = Lf[bait_cols].mean(axis=1) - Lf[ctrl_cols].mean(axis=1)
    worst = Lf[bait_cols].min(axis=1) - Lf[ctrl_cols].max(axis=1)  # weakest bait rep vs best control
    out = pd.DataFrame({'n_bait': n_bait, 'n_ctrl': n_ctrl, 'n_ctrl_runs_used': len(ctrl_cols) - len(dead),
                        'log2_enrichment': enrich, 'worst_case_log2': worst})
    out['interactor'] = (n_bait >= min_bait_reps) & ((n_ctrl == 0) | (enrich >= fc_cutoff))
    return out.sort_values('log2_enrichment', ascending=False)
```

Seeded matrix (sticky binders excluded, dead control announced, `min_bait_reps` dropout): `examples/lfq_normalization.py`.

This fold-change/presence score is the honest ceiling for ONE bait with a few controls. For a
probability rather than a cutoff use SAINTexpress (spectral counts, `interaction`/`prey`/`bait`
files -> AvgP, report BFDR <= 0.01-0.05) or CompPASS WD scores across a bait MATRIX; both need
several independent baits, or the CRAPome as an external control set, before their statistics mean
anything. Feed them counts or raw intensities -- never a median/SL/IRS-normalized matrix.
