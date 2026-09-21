# Oligo Design for Pooled Synthesis

## Oligo Design for Pooled Synthesis

**Goal:** Generate the final oligo sequence ready for chip-based synthesis. Vendor limits differ: Twist oligo pools cap at ~300 nt per oligo with no fixed pool size, GenScript's 92K format spans 20-170 nt, and Agilent OLS 244K spans 30-230 nt.

**Approach:** Add subpool PCR primers (so multiple sublibraries can share a synthesis array), the BsmBI/Esp3I overhang for golden-gate cloning into LentiGuide-Puro (Addgene 52963) or LentiCRISPRv2, and append the tracrRNA scaffold if the array length permits.

```python
def build_oligo(spacer, vector='lentiGuide-Puro', subpool_idx=None):
    '''Construct final oligo for pooled synthesis.

    LentiGuide-Puro / LentiCRISPRv2 use BsmBI (Esp3I) with these overhangs:
        forward: 5'-CACCG[spacer]-3'
        reverse: 5'-AAAC[revcomp(spacer)]C-3'
    For chip synthesis, the spacer is flanked by subpool-specific PCR primers.'''
    subpool_fwd = {
        1: 'GGAAAGGACGAAACACCG',   # subpool 1 forward primer + BsmBI overhang
        2: 'GAGGCACTGGGCAGGTACCG',
    }.get(subpool_idx, 'GGAAAGGACGAAACACCG')
    # First 33 nt of the Chen 2013 sgRNA(F+E) optimized scaffold. NOTE: lentiGuide-Puro (#52963)
    # and lentiCRISPRv2 (#52961) carry the ORIGINAL scaffold; F+E belongs to lentiCRISPRv2-Opti (#163126).
    scaffold_short = 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG'
    oligo = subpool_fwd + spacer + scaffold_short
    if len(oligo) > 200:
        raise ValueError(f'Oligo length {len(oligo)} exceeds the 200 nt design budget; check the vendor limit')
    return oligo
```

**Subpool design:** A large synthesis pool can be partitioned into multiple sublibraries via subpool primers; each sub-PCR amplifies its subpool, allowing one synthesis batch to serve several screens. Typical subpool size: 10k-20k oligos.
