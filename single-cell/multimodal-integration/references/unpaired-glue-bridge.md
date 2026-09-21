## Unpaired / Diagonal: GLUE (Python)

**Goal:** Align independent scRNA and scATAC with no shared cells via a prior feature graph.

**Approach:** Configure each dataset with a count-appropriate probabilistic model, build a gene-anchored guidance graph, fit GLUE, then read aligned embeddings.

```python
import scglue

scglue.models.configure_dataset(rna, 'NB', use_highly_variable=True, use_rep='X_pca')     # NB needs RAW counts
scglue.models.configure_dataset(atac, 'ZINB', use_highly_variable=True, use_rep='X_lsi')
graph = scglue.genomics.rna_anchored_guidance_graph(rna, atac)     # peak-near-gene prior; coords must share genome build
# GLUE also trains a VAE, so the seed is pinned explicitly for reproducibility (scglue's
# documented default is already random_seed=0; the totalVI-sized drift was not measured for
# GLUE). Checked against scglue's documented API, not run -- scglue has no Windows build.
glue = scglue.models.fit_SCGLUE({'rna': rna, 'atac': atac}, graph, init_kws={'random_seed': 0})
rna.obsm['X_glue'] = glue.encode_data('rna', rna)
atac.obsm['X_glue'] = glue.encode_data('atac', atac)
```

Verify cell-type structure is preserved (not just modality overlap); adversarial alignment can over-mix distinct populations.

## Unpaired / Diagonal: Seurat v5 Bridge Integration (R)

**Goal:** Map an unpaired scATAC query onto a labelled scRNA reference, using a paired multiome dataset as the bridge.

**Approach:** Preprocess the three datasets in their native pipelines (query ATAC only TF-IDF), build the bridge reference from the scRNA reference plus the multiome bridge, find anchors by projecting the query into the bridge's ATAC LSI space, then transfer labels and project onto the reference UMAP.

```r
library(Seurat)
library(Signac)

# rna:   labelled scRNA reference (meta.data$celltype), NormalizeData/ScaleData/RunPCA, RunUMAP(return.model = TRUE)
# multi: paired multiome bridge with 'RNA' (normalized) and 'ATAC' (RunTFIDF, RunSVD -> 'lsi') assays
# atac:  unpaired scATAC query on the SAME peak set as the bridge's ATAC assay, RunTFIDF only
bridge <- PrepareBridgeReference(
    reference = rna, bridge = multi,
    reference.reduction = 'pca', reference.dims = 1:20,
    normalization.method = 'LogNormalize',          # 'SCT' if the reference and bridge RNA were SCTransformed
    bridge.ref.assay = 'RNA', bridge.query.assay = 'ATAC',
    supervised.reduction = 'slsi', laplacian.reduction.dims = 1:20)

# dims start at 2: drop LSI_1 ONLY if DepthCor confirms it tracks depth (see `references/multiome-mofa.md`)
anchors <- FindBridgeTransferAnchors(extended.reference = bridge, query = atac,
                                     reduction = 'lsiproject', dims = 2:20)

# reference = the object PrepareBridgeReference RETURNED, not the original scRNA object
# (the anchorset lives in its 'Bridge' assay; passing the original errors "assay ... does not match")
atac <- MapQuery(anchorset = anchors, reference = bridge, query = atac,
                 refdata = list(celltype = 'celltype'), reduction.model = 'umap')
# predicted.celltype, predicted.celltype.score, and ref.umap are now on the query
```

Checked on Seurat 5.5.0 / Signac 1.17.1 with synthetic three-population data (300 reference RNA cells, 300 bridge multiome cells, 300 ATAC-only query cells, planted marker genes and peaks): 81% of query cells recovered their true type with `LogNormalize` (chance is 33%), 75% with `SCT`. Check `predicted.celltype.score` before trusting the transferred labels.
