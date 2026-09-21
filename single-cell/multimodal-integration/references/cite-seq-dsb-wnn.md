## CITE-seq: Denoise ADT, Then Joint Embed (Seurat)

**Goal:** Remove ADT background with DSB before WNN, because WNN does not denoise protein.

**Approach:** Estimate ambient from empty droplets and per-cell technical noise from a mixture plus isotype controls, then feed denoised ADT into the standard PCA -> WNN flow.

```r
library(dsb)
library(Seurat)

raw <- Read10X('raw_feature_bc_matrix/')           # unfiltered: contains empty droplets
cells <- Read10X('filtered_feature_bc_matrix/')    # called cells

adt_cells <- as.matrix(cells[['Antibody Capture']])
adt_empty <- as.matrix(raw[['Antibody Capture']][, setdiff(colnames(raw[['Antibody Capture']]), colnames(adt_cells))])

# Guard against DSB's own documented failure mode (see Common Errors): a filtered/cell
# matrix passed as empty_drop_matrix produces a plausible-looking but meaningless output
# with NO error or warning from DSBNormalizeProtein itself (verified, dsb 2.0.1). True
# empty droplets carry mostly ambient signal, so their total ADT counts must be markedly
# lower than in called cells.
med_cells <- median(colSums(adt_cells))
med_empty <- median(colSums(adt_empty))
if (med_empty >= med_cells * 0.5) {
    stop(sprintf(
        "empty_drop_matrix does not look like empty droplets (median total ADT %.1f vs cells %.1f) -- DSB needs the raw/unfiltered matrix's non-cell barcodes, not a second cell matrix.",
        med_empty, med_cells))
}

# isotype.control.name.vec must name the ACTUAL isotype rows (often IgG1/IgG2a/Mouse-IgG2b-Ctrl); the regex below misses those
# When isotypes are absent or not matched, set use.isotype.control = FALSE (keep denoise.counts = TRUE) and pass real names explicitly
adt_dsb <- DSBNormalizeProtein(
    cell_protein_matrix = adt_cells,
    empty_drop_matrix = adt_empty,
    denoise.counts = TRUE,
    use.isotype.control = TRUE,
    isotype.control.name.vec = grep('[Ii]sotype|IgG', rownames(adt_cells), value = TRUE)
)
```

## CITE-seq: WNN Joint Clustering (Seurat)

**Goal:** Build one weighted-NN graph from denoised RNA and ADT and cluster on it.

**Approach:** Reduce each modality independently (PCA on RNA, PCA on the small ADT panel), then learn per-cell modality weights and cluster/embed on the joint graph.

```r
obj[['ADT']] <- CreateAssay5Object(data = adt_dsb)        # DSB output is already normalized data
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA(reduction.name = 'pca')

DefaultAssay(obj) <- 'ADT'
VariableFeatures(obj) <- rownames(obj[['ADT']])
obj <- ScaleData(obj) |> RunPCA(reduction.name = 'apca', npcs = min(18, nrow(obj[['ADT']]) - 1))

# dims.list matched to informative dims; small ADT panels saturate by ~1:18
obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'apca'), dims.list = list(1:30, 1:18))
obj <- FindClusters(obj, graph.name = 'wsnn', algorithm = 3)   # algorithm 3 = SLM (the tutorial choice), NOT Leiden
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')

# Inspect the per-cell weight distribution; a single dominant modality is a red flag
VlnPlot(obj, features = 'RNA.weight', group.by = 'seurat_clusters')
```
