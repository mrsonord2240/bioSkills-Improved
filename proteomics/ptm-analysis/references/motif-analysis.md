# Motif Analysis with the Correct Background

## Motif Analysis with the Correct Background

**Goal:** Find kinase/writer motifs around the modified residue without rediscovering amino-acid composition bias.

**Approach:** Use the `Sequence window` (+/-15 residues, 31-mer) MaxQuant already provides, centered on the site. The background MUST be an experiment-matched S/T/Y set drawn from the identified proteins (or a central-residue-preserving shuffle), NOT the whole proteome or IUPAC-random -- those just report the composition of phospho-rich disordered regions. motif-x and MoMo p-values are only valid when the background is built this way.

```python
from collections import Counter
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

# 'Sequence window' is a 31-mer (+/-15) centered on the modified residue.
WINDOW_HALF = 7  # +/-7 flanking is the standard kinase-motif window
# foreground: e.g. the regulated-up class-I sites; here all class-I sites from the expansion block
foreground = [w[15 - WINDOW_HALF: 16 + WINDOW_HALF] for w in phospho['Sequence window'].dropna() if len(w) >= 31]

def matched_background(fasta, accessions, residues='ST'):
    '''Every S/T (or Y) window in the IDENTIFIED proteins, central residue preserved.
    fasta: {accession: sequence}; accessions: proteins identified in this experiment.'''
    windows = []
    for acc in accessions:
        seq = fasta[acc]
        for i, aa in enumerate(seq):
            if aa in residues:
                windows.append(''.join(seq[j] if 0 <= j < len(seq) else '_' for j in range(i - WINDOW_HALF, i + WINDOW_HALF + 1)))
    return windows

def position_frequencies(windows):
    counts = {i: Counter() for i in range(-WINDOW_HALF, WINDOW_HALF + 1)}
    for w in windows:
        for offset, aa in zip(range(-WINDOW_HALF, WINDOW_HALF + 1), w):
            if aa not in '_X':
                counts[offset][aa] += 1
    return counts

def motif_enrichment(fg_windows, bg_windows):
    '''One-sided Fisher test per (position, residue), BH-adjusted across all tests.'''
    fg, bg = position_frequencies(fg_windows), position_frequencies(bg_windows)
    rows = []
    for offset in fg:
        if offset == 0:
            continue
        n_fg, n_bg = sum(fg[offset].values()), sum(bg[offset].values())
        for aa, k in fg[offset].items():
            _, p = fisher_exact([[k, n_fg - k], [bg[offset][aa], n_bg - bg[offset][aa]]], alternative='greater')
            rows.append({'offset': offset, 'aa': aa, 'fg': k, 'fg_total': n_fg, 'bg': bg[offset][aa], 'bg_total': n_bg, 'p': p})
    table = pd.DataFrame(rows)
    table['q'] = false_discovery_control(table['p'], method='bh')
    return table.sort_values('p')
```

For a publication-grade enrichment logo, hand the foreground and a matched background to a dedicated tool (motif-x / MoMo) and render with data-visualization/sequence-logos.
