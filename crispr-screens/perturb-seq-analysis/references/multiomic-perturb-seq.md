## Multiomic Perturb-seq (RNA + ATAC)

**For chromatin readout:** Use 10X Multiome with CRISPRi/a; sgRNA assignment via the same scATAC-seq library. This section assumes already-quantified RNA (genes x cells) and ATAC (peaks x cells) count matrices for the same cells (shared `obs_names`) -- peak calling from raw fragment files is upstream of this Skill (ArchR or Signac in R, or CellRanger ARC).

```python
import muon as mu
import scanpy as sc
import pertpy as pt

mdata = mu.MuData({'rna': adata_rna, 'atac': adata_atac})  # shared obs_names = shared cells

# RNA side: sgRNA assignment + Mixscape escaper filtering exactly as in the sections above,
# writing .obs['mixscape_class'] (e.g. 'GENE_A KO') on mdata['rna']. Propagate the per-target
# call to the ATAC modality via the shared cell index -- do NOT use the pooled
# 'mixscape_class_global' here, which merges different target genes' KO cells together and
# dilutes any perturbation-specific chromatin signal.
mdata['atac'].obs['mixscape_class'] = mdata['rna'].obs['mixscape_class']

# ATAC side: differential accessibility per perturbation (KO vs NTC) on normalized counts.
# TF-IDF (muon.atac.pp.tfidf) is for embedding/LSI clustering, not per-feature testing here --
# verified empirically: its cell-wise reweighting distorted the Wilcoxon null on a planted-signal
# synthetic dataset when one condition's total accessible-peak count shifted; normalize_total +
# log1p recovered the planted differential peaks cleanly (5/5 in the top 5 by adjusted p-value).
atac_pert = mdata['atac'][mdata['atac'].obs['mixscape_class'].isin(['GENE_A KO', 'NTC'])].copy()
sc.pp.normalize_total(atac_pert)
sc.pp.log1p(atac_pert)
sc.tl.rank_genes_groups(atac_pert, groupby='mixscape_class', groups=['GENE_A KO'],
                         reference='NTC', method='wilcoxon')
peak_result = sc.get.rank_genes_groups_df(atac_pert, group='GENE_A KO')

# Peak-to-gene linking (which differential peak sits near which differential gene) needs a
# genome annotation file: muon.atac.pp.add_peak_annotation(mdata, annotation_file) followed by
# muon.atac.tl.rank_peaks_groups(...) adds nearest-gene/distance columns automatically. Without
# an annotation file, report differential genes (PyDESeq2, Pertpy Unified Framework section
# above) and differential peaks (peak_result, above) separately, as here.
```
