# Sequence Filtering

Read when the request is to subset sequences by ID pattern, gap content or uniqueness. Uses `normalize_alignment()` from SKILL.md (Gap and Case Normalisation).

## Sequence Filtering

**Goal:** Subset an alignment to retain only sequences matching specific criteria (ID pattern, gap content, uniqueness).

**Approach:** Iterate over alignment records, apply filter conditions, and reconstruct a new MultipleSeqAlignment from matching records. A filter that would remove every sequence raises `ValueError` instead of returning an empty alignment.

```python
import re

def _require_kept(kept, alignment, what):
    if not kept:
        raise ValueError(f'{what} removes all {len(alignment)} sequences')
    return MultipleSeqAlignment(kept, annotations=alignment.annotations,
                                column_annotations=alignment.column_annotations)

def filter_by_id(alignment, pattern):
    regex = re.compile(pattern)
    return _require_kept([r for r in alignment if regex.search(r.id)], alignment, f'pattern {pattern!r}')

def filter_by_gap_content(alignment, max_gap_fraction=0.1):
    fractions = [str(r.seq).count('-') / len(r.seq) for r in normalize_alignment(alignment)]
    kept = [r for r, f in zip(alignment, fractions) if f <= max_gap_fraction]
    return _require_kept(kept, alignment, f'max_gap_fraction={max_gap_fraction} (lowest fraction {min(fractions):.2f})')

def remove_duplicates(alignment):
    # Rows are compared after normalisation (AC-GT, AC.GT and ac-gt are one sequence); the original records are kept.
    seen, kept = set(), []
    for record, normalized in zip(alignment, normalize_alignment(alignment)):
        key = str(normalized.seq)
        if key not in seen:
            seen.add(key)
            kept.append(record)
    return _require_kept(kept, alignment, 'remove_duplicates')
```
