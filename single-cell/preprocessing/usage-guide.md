# Single-Cell Preprocessing - Usage Guide

## Overview

This skill covers quality control, ambient-RNA handling, normalization, and feature selection for single-cell RNA-seq in Scanpy (Python) and Seurat (R). It is framed around the decisions that drive every downstream result: where to set QC thresholds without deleting real cell types, whether to remove ambient RNA and with which tool, which normalization to use, and whether to scale and regress out covariates.

## Prerequisites

**Python (Scanpy):**
```bash
pip install scanpy matplotlib scikit-misc        # scikit-misc is required for flavor='seurat_v3'
```

**R (Seurat + ambient/normalization):**
```r
install.packages('Seurat')
BiocManager::install(c('scran', 'celda'))         # scran normalization; celda::decontX ambient removal
install.packages('SoupX')
# CellBender is a Python GPU CLI: pip install cellbender
```

## Quick Start

Tell the AI agent what is needed:
- "Run MAD-based QC and filter low-quality cells without deleting high-mito cell types"
- "Remove ambient RNA, then normalize and find highly variable genes"
- "Normalize with shifted-log and select 2000 HVGs from raw counts"

## Example Prompts

### Quality Control
> "Compute QC metrics with mito, ribo, and hemoglobin fractions and show the distributions"

> "Filter cells using 5 MAD on counts/genes and a tissue-aware mitochondrial cap"

> "My tissue is cardiac muscle - set a mito threshold that does not delete cardiomyocytes"

### Ambient RNA
> "Run SoupX on my Cell Ranger output and report the contamination fraction"

> "Decide whether I need CellBender or DecontX for this snRNA-seq dataset"

### Normalization and Features
> "Normalize with shifted-log and explain whether I should use scran instead"

> "Select highly variable genes from raw counts with seurat_v3"

> "Should I scale and regress out total_counts before PCA?"

## What the Agent Will Do

1. Annotate mito/ribo/hemoglobin gene sets and compute joint QC metrics
2. Set adaptive (MAD) thresholds and adjust the mito cutoff to the tissue
3. Decide whether ambient-RNA removal is warranted and pick one tool on the raw matrix
4. Stash raw counts, then normalize (shifted-log by default; scran/Pearson when justified)
5. Select HVGs with the correct input type for the chosen flavor
6. Skip or apply scaling/regression based on whether a covariate is confounded with biology
7. Run PCA on the HVG matrix, ready for clustering

## Related Skills

- single-cell/data-io - load the raw matrix before preprocessing
- single-cell/doublet-detection - per-sample doublet calling around the QC step
- single-cell/clustering - PCA, neighbors, and clustering after preprocessing
- single-cell/batch-integration - correct batch effects instead of regressing them out
- single-cell/markers-annotation - find markers after clustering
- differential-expression/deseq2-basics - pseudobulk DE across samples (avoids single-cell pseudo-replication)
