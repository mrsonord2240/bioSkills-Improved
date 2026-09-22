# Library QC and merge

Moved from SKILL.md. Read when summarising a library or merging several libraries.

## QC and Merge Libraries

**Goal:** Summarize a library and combine multiple libraries without dropping legitimate distinct transitions.

**Approach:** Report precursor/protein counts and transitions-per-precursor, then dedup on the FULL transition key. Deduping on (sequence, fragment-type, fragment-number) alone drops real transitions that differ only in precursor charge or fragment charge -- key on all five.

```python
import pandas as pd

TRANSITION_KEY = ['ModifiedSequence', 'PrecursorCharge', 'FragmentType',
                  'FragmentSeriesNumber', 'FragmentCharge']  # full key; charges matter

def merge_libraries(libs):
    combined = pd.concat(libs, ignore_index=True)
    combined['precursor_total'] = combined.groupby(
        ['ModifiedSequence', 'PrecursorCharge'])['LibraryIntensity'].transform('sum')
    combined = combined.sort_values('precursor_total', ascending=False)
    combined = combined.drop_duplicates(subset=TRANSITION_KEY).drop(columns='precursor_total')
    return combined

def library_stats(lib):
    n_prec = lib.groupby(['ModifiedSequence', 'PrecursorCharge']).ngroups
    return {'precursors': n_prec, 'proteins': lib['ProteinId'].nunique(),
            'transitions_per_precursor': round(len(lib) / n_prec, 1)}
```
