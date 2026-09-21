# Interpreting BAGEL2 calls (essential, neutral, tumor suppressor)

> Moved out of SKILL.md. "Below" and "Failure Modes" in the text refer to sections of SKILL.md.

## Interpret BAGEL2 Results

**Goal:** Stratify genes into essential, non-essential, and tumor-suppressor categories.

**Approach:** Apply BF threshold to classify essentials; only classify negative-BF genes as tumor suppressors when the screen design actually expects enrichment (see Failure Modes below) (verified below).

```python
import pandas as pd
import warnings

# Assay-control pseudo-genes spiked into CRISPR libraries (never real biology) --
# verified present and dominating the naive tumor-suppressor call on real HAP1 TKOv3
# output; extend this set to match your library's own controls.
ASSAY_CONTROLS = {'LacZ', 'luciferase', 'EGFP'}

def interpret_bagel(bf_path, bf_essential=6, bf_tumor_suppressor=-6,
                     screen_type='dropout', control_genes=ASSAY_CONTROLS,
                     tumor_suppressor_frac_warn=0.05):
    '''Classify genes from BAGEL2 BF output.

    screen_type: 'dropout' (default) only calls `essential`. Tumor-suppressor calls
    require screen_type='enrichment' or 'both' -- per the Failure Modes section below,
    a pure dropout screen's negative-BF genes are noise, not tumor suppressors.
    '''
    df = pd.read_csv(bf_path, sep='\t')
    df = df[~df['GENE'].isin(control_genes)].copy()   # drop assay-control pseudo-genes
    df['call'] = 'neutral'
    df.loc[df['BF'] > bf_essential, 'call'] = 'essential'
    if screen_type in ('enrichment', 'both'):
        df.loc[df['BF'] < bf_tumor_suppressor, 'call'] = 'tumor_suppressor'
        frac = (df['call'] == 'tumor_suppressor').mean()
        if frac > tumor_suppressor_frac_warn:
            warnings.warn(
                f"{frac:.1%} of genes flagged tumor_suppressor -- implausibly high; "
                "this usually means a dropout-only screen is being scored for "
                "enrichment. Re-check screen_type and BF<-6 calls against literature "
                "before reporting."
            )
    return df.sort_values('BF', ascending=False)
```

Verified on real HAP1 TKOv3 data (a T0-vs-T18 dropout screen): skipping the guard flags **86.5% of the genome** as tumor-suppressor, with three assay-control pseudo-genes (`LacZ`, `luciferase`, `EGFP`) as the top hits. With the guard, the default
(`screen_type='dropout'`) returns 0 tumor-suppressor calls and excludes the 3 control
pseudo-genes; explicit `screen_type='enrichment'` still surfaces the true tumor
suppressors TSC1/TSC2 as the two most-negative-BF entries (now that the control genes
that previously masked them are excluded) but raises the 86.6%-flagged warning so the
caller doesn't report it uncritically.

**Tumor suppressor identification:** Genes with significantly negative BF (e.g., <-6) in a screen actually designed to detect enrichment (drug-resistance, GoF) indicate fitness advantage from their loss, which is biologically distinct from "non-essential". **Do not call tumor suppressors from a pure dropout screen** -- see Failure Modes.
