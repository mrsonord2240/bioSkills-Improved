# CRISPRi / CRISPRa TSS Targeting

## CRISPRi / CRISPRa TSS Targeting

**Goal:** Position guides relative to the empirical TSS for maximum knockdown (CRISPRi) or activation (CRISPRa).

**Approach:** Resolve TSS from FANTOM5 CAGE peaks (highest-ranked peak per gene; fall back to Ensembl/RefSeq if absent), define the modality-specific window, score candidate spacers in that window with Rule Set 2 plus the Horlbeck/Sanson CRISPRi/a-tailored rules, and select 5-6 guides per gene biased toward the window center.

```python
def crispri_window(tss_coord, strand='+'):
    '''Dolcetto convention: search -50 to +300 around the FANTOM5 highest-rank CAGE peak.
    Reason: Sanson 2018 found +25 to +75 nt downstream of the TSS optimal for CRISPRi,
    so rank candidates toward that band; the search is relaxed outward to fill the
    per-gene guide quota when poorly-annotated TSSs leave too few candidates.'''
    if strand == '+':
        return (tss_coord - 50, tss_coord + 300)
    return (tss_coord - 300, tss_coord + 50)

def crispra_window(tss_coord, strand='+'):
    '''Calabrese convention: -150 to -75 upstream of TSS.
    Reason: dCas9-VP64 (and SAM, SunTag) activate maximally when bound
    just upstream of Pol II loading. Horlbeck v2 CRISPRa uses -550 to -25
    (broader, lower per-guide signal). For SAM, prefer Calabrese tightness;
    for SunTag, Horlbeck width is acceptable.
    Caveat: at only 75bp wide, this window routinely fails to contain a full
    6-guide quota's worth of PAM sites passing the GC/poly-T filter -- budget
    for shortfalls (report actual count per gene rather than padding with
    out-of-window guides) or widen to Horlbeck v2 when the quota must be met.'''
    if strand == '+':
        return (tss_coord - 150, tss_coord - 75)
    return (tss_coord + 75, tss_coord + 150)
```

**Critical nuance:** Cell-type-specific TSSs differ from the FANTOM5 consensus in ~15% of genes. For tissue-specific screens (e.g., neuron, hepatocyte), re-derive TSSs from a matched CAGE / GRO-seq / PRO-seq dataset before locking guide positions, or knockdown efficiency drops several-fold. The single most common cause of "weak" CRISPRi hits is mis-positioned guides against an alternative TSS: dCas9-KRAB knockdown is maximal within +/-100 bp of the actual Pol II loading site, and canonical Ensembl/RefSeq annotation can be off by 1-10 kb for genes with broad or non-canonical promoters. Symptom: "easy" essentials (RPS, RPL, EIF) show normal dropout but newer genes do not, and the library validates poorly against CEGv2.
