# LCV (Latent Causal Variable)

## LCV (Latent Causal Variable)

LCV uses LDSC-merged genome-wide sumstats and reports gcp (genetic causality proportion) on [-1, 1]. It is a complement to, not a replacement for, MR; gcp ~ 0 with high LDSC rg implies pure genetic correlation without partial causation. LCV uses ALL genome-wide SNPs after LDSC-merging, not the MR instrument set.

```r
source('LCV/R/RunLCV.R')
res_lcv <- RunLCV(ldscores$L2, x$Z, y$Z)
# res_lcv$gcp.pm (posterior mean gcp; there is no res_lcv$gcp field); res_lcv$pval.gcpzero.2tailed
```

**gcp interpretation:**

| gcp value | Interpretation |
|-----------|----------------|
| 0 | Pure genetic correlation; no partial causation |
| 0.5 | Partial causation; mixture |
| 0.6 | Partial causation; modestly causal direction |
| 1 | Fully causal in tested direction |
| Significant p_gcp != 0 | Directional evidence of (partial) causation |
