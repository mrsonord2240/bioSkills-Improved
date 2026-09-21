# DIA-NN workflow

Read this when the search output is DIA-NN `report.parquet`: q-value filtering before pivoting, the PG.MaxLFQ matrix, and the 0 -> NA step before log2.

### DIA-NN Workflow
DIA-NN 1.9+ defaults to report.parquet (the only default in 2.0); read it with arrow, not read.delim. Filter on q-values BEFORE pivoting, or low-confidence rows enter the matrix. Route to proteomics/dia-analysis for the mechanics.
```r
library(arrow)
library(dplyr)
library(tidyr)

diann <- read_parquet('report.parquet')

# Filter to 1% FDR at precursor AND protein-group level before pivoting.
# Use the GLOBAL protein-group q-value for the cross-run matrix (per-run min(Q.Value) is anti-conservative).
# When MBR is ON, MBR has its own FDR: add the Lib.* q-values (Lib.Q.Value, Lib.PG.Q.Value <= 0.01).
diann_filt <- diann %>%
    filter(Q.Value <= 0.01 & PG.Q.Value <= 0.01 & Global.PG.Q.Value <= 0.01)

# PG.MaxLFQ is ALREADY cross-run MaxLFQ-normalized at report generation. Re-normalizing it
# double-normalizes -- go straight to log2 + limma with no further normalization. To apply the
# skill's own median-centering instead, pivot raw PG.Quantity here, not PG.MaxLFQ.
protein_matrix <- diann_filt %>%
    select(Protein.Group, Run, PG.MaxLFQ) %>%
    distinct() %>%
    pivot_wider(names_from = Run, values_from = PG.MaxLFQ)

# PG.MaxLFQ path: log2-transform and go straight to limma (no re-normalization). DIA-NN writes 0
# for "not quantified in this run", so 0 -> NA FIRST: log2(0) is -Inf, and -Inf cells propagate
# silently until eBayes stops with "missing value where TRUE/FALSE needed".
m <- as.matrix(protein_matrix[, -1])
rownames(m) <- protein_matrix$Protein.Group
m[m == 0] <- NA
log2_matrix <- log2(m)
stopifnot(!any(is.infinite(log2_matrix)))
```
