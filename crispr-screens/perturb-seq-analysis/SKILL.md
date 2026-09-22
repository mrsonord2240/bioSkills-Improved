---
name: bio-crispr-screens-perturb-seq-analysis
description: Analyzes single-cell pooled CRISPR screens (Perturb-seq, CROP-seq, Perturb-CITE-seq, ECCITE-seq, multiome) where each cell carries an sgRNA and a scRNA-seq / surface-protein / chromatin readout. Covers experimental design (direct-capture Perturb-seq Dixit 2016 vs CROP-seq 3'UTR-barcoded Datlinger 2017 vs ECCITE-seq vs Multiome), MOI for sgRNA assignment, escaper-cell filtering (Mixscape, Papalexi 2021), SCEPTRE NB GLM + permutation for low-MOI (Barry 2024 Genome Biol 25:124), the Pertpy framework, factor decomposition, genome-scale Perturb-seq (Replogle 2022 Cell, 2.5M cells), and per-perturbation single-cell DE. Use when running a single-cell CRISPR screen, choosing direct-capture vs CROP-seq architecture, filtering escaper cells, performing single-cell DE, integrating Perturb-seq with pathway analysis, scaling to GW CRISPRi via Replogle protocol, or analyzing multi-omics screens.
tool_type: python
primary_tool: Pertpy
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples checked on pertpy 1.3.0, scanpy 1.12+, anndata 0.13+, pandas 2.2+, numpy 1.26+, scipy 1.12+, sceptre 0.99.0 (R / katsevich-lab/sceptre GitHub HEAD), muon 0.1.9.

**pertpy >= 1.0 breaking changes from the 0.6.x examples some agents may have seen:** `PyDESeq2.test_contrasts()` takes a numeric contrast vector built via `de.contrast(column, baseline, group_to_compare)`, not a `contrast=(column, group, baseline)` tuple; and result columns are `log_fc` / `p_value` / `adj_p_value`, not `log2FoldChange` / `padj`. All code blocks below already use the current API.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show pertpy scanpy anndata`
- R: `packageVersion('sceptre')`; `?sceptre`; `?Seurat::PrepLDA`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Single-Cell Perturb-Seq Analysis

**"Analyze a single-cell pooled CRISPR perturbation screen"** -> Assign sgRNAs to cells, filter unperturbed escapers, normalize counts, fit per-gene differential expression conditioned on perturbation, and rank perturbations by their molecular effect.

- Python: `pertpy` unified framework for Mixscape + SCEPTRE-via-R + differential expression
- R: `sceptre` for low-MOI NB GLM + permutation testing
- Python/R: `Seurat::MixscapeLDA` and downstream

## Experimental Architecture Comparison

| Method | Year | Architecture | Readout | MOI | Single-cell sgRNA detection |
|--------|------|--------------|---------|-----|------------------------------|
| Perturb-seq (Dixit 2016, *Cell*) | 2016 | sgRNA expressed in cassette; direct PCR capture | scRNA-seq | Low-to-moderate (MOI ~0.35-1.4; a minority of cells receive multiple guides, enabling epistasis analysis) | Yes via amplicon-PCR pre-sequencing |
| CROP-seq (Datlinger 2017, *Nat Methods*) | 2017 | hU6-sgRNA cassette placed in the 3' LTR of lentiGuide-Puro; LTR duplication puts the sgRNA in the 3'UTR of the Pol II puromycin-resistance transcript | scRNA-seq | Low (1-2 sgRNAs/cell) | Native via 10X 3' chemistry |
| Perturb-CITE-seq (Frangieh 2021, *Nat Genet*) | 2021 | Adds surface-protein hashtag oligos to CROP-seq | scRNA-seq + ADT (protein) | Low | CROP-seq architecture |
| ECCITE-seq (Mimitou 2019, *Nat Methods*) | 2019 | Surface-protein hashtag with sgRNA-marked cells | scRNA-seq + ADT | Low | Hash + sgRNA |
| Perturb-ATAC (Rubin 2019, *Cell*) | 2019 | scATAC-seq readout | scATAC | Low | sgRNA capture via separate library prep |
| Perturb-multiome (10X) | 2021+ | scRNA + scATAC simultaneously | scRNA + ATAC | Low | Direct capture from sgRNA cassette |
| Replogle GW Perturb-seq (2022, *Cell*) | 2022 | Multiplexed CRISPRi with sgRNA barcoding | scRNA-seq | 1 sgRNA/cell | Direct capture |

**Decision rule:** Standard scRNA + sgRNA at low cost -> CROP-seq. Genome-wide CRISPRi screens -> Replogle's CRISPRi + 10X 3' direct-capture protocol, the gold standard (>2.5M cells; the genome-scale K562 screen targeted ~9,866 expressed genes in Replogle 2022). Protein readout -> Perturb-CITE-seq. Chromatin readout -> Perturb-multiome. Hashed cells + sgRNA -> ECCITE-seq. Low-throughput pilot -> original Dixit Perturb-seq. Genome-scale design and budget: `references/genome-wide-perturb-seq.md`; chromatin readout analysis: `references/multiomic-perturb-seq.md`.

## MOI and sgRNA Assignment

**The central technical challenge:** Each cell must receive exactly one sgRNA (otherwise the perturbation is undefined). At MOI 0.3 (the standard target), ~26% of cells get ≥1 sgRNA, but 4% get ≥2; at MOI 0.5, ~9% of cells get multiple sgRNAs. The cells with multiple sgRNAs must be filtered or analyzed as combinatorial perturbations (see crispr-screens/combinatorial-screens for intentionally-high-MOI paired-guide design).

**Assignment workflow:**

1. **Detect sgRNA reads per cell:** From the sgRNA library prep (direct capture or 3'UTR barcode), count reads per sgRNA per cell.
2. **Threshold:** Most pipelines use 10+ reads of one sgRNA to assign that perturbation.
3. **Multiplets:** Cells with 2+ sgRNAs at >10 reads each are either multi-perturbed (analyzable as combinatorial) or doublets.
4. **Doublet detection:** Use scDblFinder, Scrublet, or AMULET (multiome) to identify doublets independently from sgRNA assignment.

**Goal:** Assign a single perturbation identity (or 'multiplet'/'none') to every cell from the sgRNA counts matrix.

**Approach:** Threshold per-cell sgRNA reads at ≥10 (Pertpy convention); cells exceeding the threshold for exactly one sgRNA are assigned that perturbation; cells with multiple sgRNAs above threshold are flagged as multiplets for filtering or combinatorial analysis.

```python
# sgRNA assignment via threshold counting
def assign_sgrna(adata, sgrna_counts_layer='sgrna_counts', threshold=10):
    '''Per-cell sgRNA assignment. Returns single assignment or 'multiplet'/'none'.'''
    import numpy as np
    counts = adata.layers[sgrna_counts_layer]  # cells x sgRNAs
    above_thresh = counts >= threshold
    n_sgrna_per_cell = above_thresh.sum(axis=1)
    assignments = np.where(
        n_sgrna_per_cell == 0, 'none',
        np.where(n_sgrna_per_cell == 1,
                  [adata.var_names[i] for i in counts.argmax(axis=1)],
                  'multiplet'))
    adata.obs['sgrna_assignment'] = assignments
    return adata
```

## Escaper Cell Filtering (Mixscape)

**Why this matters:** Not all sgRNA-positive cells actually edit. The escaper fraction is guide- and gene-dependent: Papalexi 2021 measured ~25% escapers for IFNGR2, perturbation rates of 39-92% across four IRF1 guides (i.e. 8-61% escapers), and no detectable perturbation at all for 15 genes. Including escapers dilutes the perturbation effect; Mixscape identifies and filters them.

**Mixscape algorithm:** For each perturbed cell, compute a "perturbation signature" = (its expression) - (mean of K nearest non-targeting-control cells). This signature isolates the perturbation effect from cell-state variation. Cells with perturbation signature similar to NTC distribution are escapers.

```python
import pertpy as pt
import scanpy as sc

# adata is a scRNA-seq AnnData with 'sgrna_assignment' column
# Pertpy 0.6+ Mixscape API (verify against installed pertpy with help(pt.tl.Mixscape))
mixscape = pt.tl.Mixscape()
mixscape.perturbation_signature(
    adata=adata,
    pert_key='sgrna_assignment',     # .obs column with sgRNA target per cell
    control='NTC',                   # match this to whatever your data's NTC label actually is
    n_neighbors=20,                  # K neighbors for KNN-NTC subtraction
    random_state=0,                  # forwarded to pynndescent.NNDescent -- omitting it makes
                                      # X_pert non-deterministic across reruns (verified: 5/2000
                                      # cells drifted between two unseeded runs on identical input)
)
# Writes .layers['X_pert'] with perturbation-signature-corrected expression

# Filter escapers: classify perturbed cells as KO (true perturbation) or NP (non-perturbed/escaper)
mixscape.mixscape(
    adata=adata,
    pert_key='sgrna_assignment',     # pert_key (not 'labels' in modern pertpy)
    control='NTC',
    new_class_name='mixscape_class', # .obs column to write
)
# Defaults to layer='X_pert' (output of perturbation_signature)

# Keep only KO cells for downstream analysis
# (mixscape_class holds '<gene> KO'; the bare label is in mixscape_class_global)
adata_ko = adata[adata.obs['mixscape_class_global'].isin(['KO'])].copy()
print(f'KO cells: {adata_ko.n_obs} ({adata_ko.n_obs/(adata.obs["sgrna_assignment"] != "NTC").sum():.1%} of perturbed)')
```

**Critical:** Mixscape can fail when the perturbation has weak phenotype; empirically Mixscape detects perturbations with log-fold-change <-0.5 (depletion) reliably, but weaker effects collapse into the NTC distribution. For genome-wide screens, run Mixscape per perturbation; for low-effect perturbations, trust the assignment without filtering.

## SCEPTRE for Low-MOI Differential Expression

SCEPTRE (NB GLM + conditional resampling, calibrated FDR) is the low-MOI DE method of choice, run in R. Full method and code: `references/sceptre-low-moi.md`.

## Pertpy Unified Framework

**Pertpy** (https://pertpy.readthedocs.io) integrates Mixscape, distance-based perturbation comparison, EdgeR/PyDESeq2/WilcoxonTest DE, and factor models in a single AnnData-based interface. For SCEPTRE specifically, invoke the R sceptre package separately (Pertpy does not wrap it).

```python
import pertpy as pt
import scanpy as sc

# Load data
mdata = pt.dt.papalexi_2021()      # returns a MuData object (rna/adt/hto/gdo modalities)
adata = mdata['rna']  # built-in example from Mixscape paper
# The per-cell target-gene label ('gene_target', with control cells labeled 'NT') lives on the
# top-level MuData.obs, not on the rna modality's own .obs -- merge it in (verified: same
# obs_names/order across mdata and mdata['rna']).
adata.obs['gene_target'] = mdata.obs['gene_target']
adata.layers['counts'] = adata.X.copy()  # PyDESeq2 needs raw counts -- save before normalizing

# Standard scRNA-seq preprocessing (scanpy)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Mixscape escaper filtering (writes .layers['X_pert'] and .obs['mixscape_class'])
# NOTE: papalexi_2021()'s own non-targeting-control label is literally 'NT' (not the 'NTC'
# placeholder used elsewhere in this Skill) -- always match `control=` to your data's actual label.
ms = pt.tl.Mixscape()
ms.perturbation_signature(adata, pert_key='gene_target', control='NT', n_neighbors=20, random_state=0)
ms.mixscape(adata, pert_key='gene_target', control='NT')

# Filter to KO cells
adata_ko = adata[adata.obs['mixscape_class_global'].isin(['KO', 'NT'])].copy()

# Pseudobulk differential expression via pertpy (PyDESeq2 backend)
# pertpy >= 1.0: build the contrast with .contrast(column, baseline, group_to_compare), then pass
# the resulting vector to test_contrasts(); result columns are log_fc / p_value / adj_p_value
# (checked on pertpy 1.3.0 -- see Version Compatibility).
# layer='counts': PyDESeq2 requires raw (near-integer) counts -- pointing it at log-normalized
# data raises "ValueError: Non-zero elements of the matrix must be close to integer values."
# (verified: reproduced this exact error on real data, then fixed it by adding the counts layer).
de = pt.tl.PyDESeq2(adata_ko, design='~gene_target', layer='counts')
de.fit()
results_df = de.test_contrasts(de.contrast('gene_target', 'NT', 'GENE_X'))

# For calibrated SCEPTRE on low-MOI single-cell data, use R sceptre directly
# (Barry 2024 Genome Biol; not bundled in pertpy)
```

## Genome-Wide Perturb-Seq (Replogle 2022)

Genome-scale CRISPRi design, cell/channel budget and the scaling formula for a different gene count: `references/genome-wide-perturb-seq.md`.

## Factor-Based Analysis

Shared-factor decomposition of perturbation effects (FR-Perturb, standalone CLI): `references/factor-decomposition.md`.

## Multiomic Perturb-seq (RNA + ATAC)

RNA + ATAC Perturb-seq (10X Multiome): propagate the Mixscape call to the ATAC modality and test differential accessibility: `references/multiomic-perturb-seq.md`.

## Reference Files

| File | Read when |
|------|-----------|
| `references/sceptre-low-moi.md` | Calibrated single-cell DE for a low-MOI screen (R `sceptre` pipeline) |
| `references/genome-wide-perturb-seq.md` | Designing or budgeting a genome-scale CRISPRi Perturb-seq (Replogle 2022) |
| `references/factor-decomposition.md` | Decomposing perturbation effects into shared factors (FR-Perturb) |
| `references/multiomic-perturb-seq.md` | RNA + ATAC (10X Multiome) Perturb-seq: differential accessibility per perturbation |

## Failure Modes

### Low sgRNA detection per cell

**Trigger:** Direct-capture method on CROP-seq library, or 3'UTR barcoding on direct-capture library.
**Mechanism:** Architecture mismatch -- the sgRNA can't be detected by the wrong library prep.
**Symptom:** sgRNA assignment rate <50% of cells.
**Fix:** Match library prep to architecture; for CROP-seq, use 10X 3' chemistry; for direct-capture Perturb-seq, use the Dixit amplicon-PCR pre-sequencing.

### Mixscape filters too many cells as escapers

**Trigger:** Weak perturbation phenotype; Mixscape's NTC-subtracted signature is similar to NTC null.
**Mechanism:** Mixscape assumes a detectable signal; weak knockdown is misclassified as escaper.
**Symptom:** >50% of perturbed cells classified as "NP" (non-perturbed); known essentials show no effect.
**Fix:** Lower Mixscape stringency; skip Mixscape for low-effect perturbations; verify Cas9 expression first.

### Doublet contamination drives apparent multi-perturbation cells

**Trigger:** High cell density loading on 10X channels.
**Mechanism:** Two cells in one droplet appear to carry two sgRNAs.
**Symptom:** "Multiplet" rate >5% after sgRNA assignment.
**Fix:** Reduce cell loading per channel (5,000-7,000 instead of 10,000); run Scrublet or scDblFinder; remove doublets before sgRNA assignment.

### MAST or Wilcoxon over-call hits

**Trigger:** Using parametric DE tools on sparse, zero-inflated scRNA-seq.
**Mechanism:** These tools assume Gaussian or simpler null; single-cell data has zero-inflation that makes them over-confident.
**Symptom:** Thousands of significant DE genes per perturbation; FDR uncalibrated.
**Fix:** Use SCEPTRE (permutation-based NB GLM); Barry 2024 benchmark shows this is the only method with calibrated FDR.

### Genome-scale Perturb-seq with insufficient cells per perturbation

**Trigger:** <500 cells per perturbation in genome-scale experiment.
**Mechanism:** DE estimation requires sufficient cells per condition; <500 lacks power for moderate effects.
**Symptom:** Inconsistent hit calls across replicates; pathway analysis non-specific.
**Fix:** Scale up cell numbers; or run focused (sub-genome) Perturb-seq with more cells per pert.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| MOI for single sgRNA per cell | 0.3 | Poisson math; ~26% infected, 4% multi-infected |
| sgRNA assignment threshold | ≥10 reads of one sgRNA | Pertpy / direct-capture convention |
| Multiplet rate (post-doublet filter) | <5% | Typical 10X 3' chemistry |
| Mixscape KO retention | Guide-dependent; 39-92% observed | Papalexi 2021 |
| Cells per perturbation (DE power) | 500-1,000 minimum genome-scale; 1,000-2,000 focused (specific module); 5,000+ single-pert deep; 2,000+ per pair combinatorial | Power convention (Replogle 2022 screened at a median >100) |
| SCEPTRE permutations | 1,000+ | Barry 2024 |
| Genes per cell (QC) | ≥500-1,000 | Standard scRNA QC |
| Mt% threshold | <15-20% | Standard scRNA QC |
| Doublet detection threshold | scDblFinder, Scrublet defaults | Methods agree |
| NTC (non-targeting control) representation | ~5% of library | Standard pooled-screen library design |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Low sgRNA detection | Architecture mismatch | Match library prep |
| Too many escapers in Mixscape | Weak phenotype | Skip Mixscape; verify Cas9 |
| Inflated DE hits | MAST / Wilcoxon used | Switch to SCEPTRE |
| Inconsistent gene effects between channels | Channel batch effect | Add channel as covariate in SCEPTRE |
| Multiplet rate >10% | Over-loading cells | Reduce loading; doublet filter |
| Per-pert DE with <100 cells | Insufficient power | Increase cell numbers; or accept low resolution |

## Scope

This Skill analyzes single-cell pooled CRISPR screen data for research purposes. A perturbation's molecular effect here (a DE gene, an escaper-filtered phenotype, a factor loading) is a research finding, not a validated therapeutic target or a patient-treatment recommendation -- a screen-level knockdown effect does not by itself establish clinical benefit or harm for a patient. Do not use this Skill's output to make or imply a diagnostic or treatment decision for an individual.

## References

- Dixit A et al. 2016. *Cell* 167:1853. Original Perturb-seq.
- Datlinger P et al. 2017. *Nat Methods* 14:297. CROP-seq.
- Frangieh CJ et al. 2021. *Nat Genet* 53:332. Perturb-CITE-seq.
- Mimitou EP et al. 2019. *Nat Methods* 16:409. ECCITE-seq.
- Rubin AJ et al. 2019. *Cell* 176:361. Perturb-ATAC.
- Papalexi E et al. 2021. *Nat Genet* 53:322. Mixscape.
- Barry T, Mason K, Roeder K, Katsevich E. 2024. *Genome Biol* 25:124. SCEPTRE for low-MOI Perturb-seq.
- Replogle JM et al. 2022. *Cell* 185:2559. Genome-wide Perturb-seq.
- Heumos L et al. 2026. *Nat Methods* 23:350-359. DOI 10.1038/s41592-025-02909-7. Pertpy framework.
- Jiang L et al. 2025. *Nat Cell Biol* 27:505. Mixscale (perturbation-strength-aware Perturb-seq).

## Related Skills

- crispr-screens/library-design - Direct-capture vs CROP-seq library design
- crispr-screens/screen-qc - sgRNA assignment rates as QC
- crispr-screens/mageck-analysis - Pseudobulk analysis as alternative
- crispr-screens/hit-calling - Pseudo-bulk hit calling alternative
- single-cell/preprocessing - scRNA-seq preprocessing
- single-cell/clustering - Post-DE clustering
- single-cell/multimodal-integration - Multiome Perturb-seq
- single-cell/perturb-seq - General single-cell screen analysis
- pathway-analysis/go-enrichment - Pathway enrichment of perturbation hits
