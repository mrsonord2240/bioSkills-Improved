# FUSION Pipeline

**Goal:** Run sumstat TWAS using FUSION and conditional joint analysis to identify independently associated genes.

**Approach:** Format GWAS sumstats to FUSION expected columns (SNP A1 A2 Z); run `FUSION.assoc_test.R` per chromosome with pre-computed weights and ancestry-matched LD; post-process with `FUSION.post_process.R` (conditional/joint analysis run by default) to identify independent genes; flag conditional-significant genes for follow-up.

```bash
# Pre-computed FUSION weights live at http://gusevlab.org/projects/fusion/
# Example: GTEx v8 Whole_Blood; download .pos summary + per-gene RData files into wgt_dir/

# GWAS sumstats expected columns: SNP A1 A2 Z (Z-score on standardised scale)
# Use TwoSampleMR or a custom munger to harmonise alleles upstream

for chr in {1..22}; do
    Rscript FUSION.assoc_test.R \
        --sumstats gwas.sumstats \
        --weights gtex_whole_blood.pos \
        --weights_dir gtex_whole_blood_wgt/ \
        --ref_ld_chr 1000G_EUR_LD/EUR. \
        --chr ${chr} \
        --out twas_chr${chr}.dat
done
cat twas_chr*.dat > twas_all.dat

# Conditional joint analysis at each significant locus
Rscript FUSION.post_process.R \
    --sumstats gwas.sumstats \
    --input twas_all.dat \
    --out twas_joint.dat \
    --ref_ld_chr 1000G_EUR_LD/EUR. \
    --chr 22 \
    --plot --locus_win 100000
# twas_joint.dat reports per-gene conditional Z; genes with joint Z > 4 are independent
```

`FUSION.post_process.R` returns conditionally independent gene signals, not PIPs; use FOCUS for probabilistic fine-mapping. The `--locus_win 100000` parameter defines the conditioning window; 100 kb is conservative for non-HLA loci. FUSION's `--coloc_P` flag runs single-SNP coloc internally but is less robust than running coloc separately on the per-gene top eQTL.

**Known upstream crash on single-SNP ("top1") genes** -- see Common Errors in SKILL.md for the exact error, cause, and a verified two-line patch.
