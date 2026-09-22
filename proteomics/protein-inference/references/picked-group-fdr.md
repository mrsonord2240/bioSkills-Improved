# Picked Protein-Group FDR

### Picked Protein-Group FDR

**Goal:** Estimate protein-group FDR without the inflation that the reused PSM formula causes on large data.

**Approach:** For each target group, find its decoy counterpart (same accessions with the decoy prefix); keep only the higher-scoring member of each target/decoy PAIR; rank the picked set and count decoys as the FDR estimate. For idXML, prefer the built-in `FalseDiscoveryRate().applyPickedProteinFDR(prot_id, String(prefix), True, True)` shown in `pyopenms-basic-inference.md`; for group-level work on large data, the kusterlab `picked_group_fdr` package implements The 2022. The sketch below pairs by the exact accession set, so decoy groups whose membership differs from their target's stay unpaired and are counted unpicked. The decoy prefix is tool-specific (`DECOY_` OpenMS/Comet, `rev_` FragPipe/Philosopher's `--tag` default, `REV__` MaxQuant) and must be passed explicitly (`DECOY_PREFIX` in the blocks). A wrong prefix does not pass silently in either implementation: `applyPickedProteinFDR` raises `IndexError: invalid unordered_map<K, T> key`, and the sketch below raises `ValueError: no decoy groups with prefix ...`. Treat both as "the prefix is wrong", not as a corrupt input file.

```python
def picked_group_fdr(groups, decoy_prefix, min_decoys=10):
    # groups: list of dicts with 'accessions' (str), 'score' (higher = better), 'is_decoy'
    n_decoy = sum(g['is_decoy'] for g in groups)
    if n_decoy == 0 or not any(a.startswith(decoy_prefix) for g in groups for a in g['accessions']):
        raise ValueError(f'no decoy groups with prefix {decoy_prefix!r}; keep decoys upstream or fix the prefix')
    if n_decoy < min_decoys:
        print(f'WARNING: only {n_decoy} decoy groups; the protein FDR estimate is not meaningful')
    by_base = {}
    for g in groups:
        base = frozenset(a.replace(decoy_prefix, '') for a in g['accessions'])
        # keep only the higher-scoring of the target/decoy pair (the 'pick')
        if base not in by_base or g['score'] > by_base[base]['score']:
            by_base[base] = g
    picked = sorted(by_base.values(), key=lambda g: g['score'], reverse=True)

    targets = decoys = 0
    for g in picked:
        if g['is_decoy']:
            decoys += 1
        else:
            targets += 1
        g['fdr'] = decoys / targets if targets else 1.0
    running_min = 1.0
    for g in reversed(picked):  # monotone q-values from the bottom up
        running_min = min(running_min, g['fdr'])
        g['qvalue'] = running_min
    return [g for g in picked if not g['is_decoy'] and g['qvalue'] <= 0.01]
```
