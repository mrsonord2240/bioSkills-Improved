# SuSiEx Cross-Ancestry Fine-Mapping

## Cross-Ancestry Fine-Mapping with SuSiEx

**Goal:** Jointly fine-map a locus across multiple ancestries assuming shared causal variants but population-specific LD.

**Approach:** Per-ancestry summary statistics + per-ancestry LD reference; SuSiEx runs a joint SuSiE model with population-specific R matrices. SuSiEx assigns populations by the ORDER of the comma-separated `--sst_file`/`--n_gwas`/`--ref_file`/`--ld_file` lists (there is no `--pop` flag); keep all four lists in the same population order.

```bash
SuSiEx \
    --sst_file=eur_sumstats.txt,eas_sumstats.txt,afr_sumstats.txt \
    --n_gwas=500000,200000,80000 \
    --ref_file=1000G_EUR,1000G_EAS,1000G_AFR \
    --ld_file=eur_ld,eas_ld,afr_ld \
    --out_dir=susiex_out \
    --out_name=locus1 \
    --chr=6 --bp=30000000,31000000 \
    --chr_col=1,1,1 --snp_col=2,2,2 --bp_col=3,3,3 \
    --a1_col=4,4,4 --a2_col=5,5,5 --eff_col=6,6,6 \
    --se_col=7,7,7 --pval_col=8,8,8 \
    --level=0.95
```

The output includes per-population PIPs and a joint credible set. Credible sets from SuSiEx are typically 2-5x smaller than EUR-only susie_rss when AFR is included, because AFR shorter LD blocks resolve EUR-tagged regions.
