# Consensus Sequence

Read when the request is for a consensus sequence. Uses `normalize_alignment()` from SKILL.md (Gap and Case Normalisation).

## Consensus Sequence

**"Get consensus sequence"** -> Derive a single representative sequence from an MSA based on majority-rule voting at each column.

**Goal:** Generate a consensus sequence from the alignment using a frequency threshold.

**Approach:** At each column, select the most common non-gap character if its share of all rows (gap rows included) reaches the threshold; otherwise mark as ambiguous. The placeholder is alphabet-aware: `N` for nucleotide, `X` for protein (`N` is asparagine and cannot be told apart from a real Asn column).

### Simple Majority Consensus
```python
def is_nucleotide(alignment, min_fraction=0.9):
    text = ''.join(str(r.seq) for r in alignment).upper().replace('-', '').replace('.', '')
    return bool(text) and sum(text.count(c) for c in 'ACGTUN') / len(text) >= min_fraction

def consensus_sequence(alignment, threshold=0.5, gap_char='-', ambiguous=None, weights=None):
    alignment = normalize_alignment(alignment)
    if ambiguous is None:
        ambiguous = 'N' if is_nucleotide(alignment) else 'X'
    weights = np.ones(len(alignment)) if weights is None else np.asarray(weights, dtype=float)
    if len(weights) != len(alignment):
        raise ValueError('weights must have one value per sequence')
    if not weights.sum() > 0:
        raise ValueError('weights must sum to a positive value')
    consensus = []
    for col_idx in range(alignment.get_alignment_length()):
        counts = Counter()
        for char, weight in zip(alignment[:, col_idx], weights):
            counts[char] += weight
        counts.pop('-', None)
        if not counts:
            consensus.append(gap_char)
            continue
        most_common_char, most_common_count = counts.most_common(1)[0]
        consensus.append(most_common_char if most_common_count / weights.sum() >= threshold else ambiguous)
    return ''.join(consensus)

consensus = consensus_sequence(alignment, threshold=0.5)
```

### Note on Bio.Align.AlignInfo
`AlignInfo.SummaryInfo` keeps only `get_column` in Biopython 1.88 (`dumb_consensus`, `gap_consensus`, `pos_specific_score_matrix`, `information_content` were removed and raise `AttributeError`). Use the custom `consensus_sequence()` above.
