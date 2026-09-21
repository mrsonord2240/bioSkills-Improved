# LDAK SumHer (reference for heritability-partitioning SKILL.md)

## LDAK SumHer Pipeline

**Goal:** Alternative h2 and functional enrichment estimate using the LDAK-Thin model for reconciliation with LDSC.

**Approach:** Reformat sumstats to LDAK input -> run `ldak --sum-hers` against the pre-computed LDAK-Thin tagging file -> compare to LDSC.

```bash
# 1. LDAK requires header: Predictor A1 A2 n Z (Z optional; can use beta + se instead)
# Reference LDAK-Thin tagging files at dougspeed.com/pre-computed-tagging-files
ldak --sum-hers trait_sumher \
    --summary trait_ldak.txt \
    --tagfile ldak.thin.hapmap.gbr.tagging \
    --check-sums NO

# 2. Partitioned with BaselineLD annotations (two steps)
# BaselineLD provides binary + continuous annotations covering coding/conserved/regulatory/MAF
# bins; download BaselineLD.zip (96 annotations) from dougspeed.com/resources and extract to ./BaselineLD/BaselineLD{1..96} (the run uses the first 86).
# The annotation flags belong to --calc-tagging (which builds the tagging file), NOT to
# --sum-hers. LDAK uses --annotation-number + --annotation-prefix (continuous) or
# --partition-number + --partition-prefix (binary). No --category-file flag exists.

# 2a. Build the annotated tagging file from a genotype reference (--power -.25 = LDAK-Thin)
ldak --calc-tagging trait_bld --bfile ref_panel --power -.25 \
    --annotation-number 86 \
    --annotation-prefix BaselineLD/BaselineLD

# 2b. Estimate partitioned h2 against that tagging file (no annotation flags here)
ldak --sum-hers trait_bld --summary trait_ldak.txt \
    --tagfile trait_bld.tagging \
    --check-sums NO
```

Pre-computed tagging files exist for GBR (HapMap reference); other ancestries require building tagging file via `--calc-tagging`. LDAK SumHer outputs h2 estimate, per-category h2 share, and enrichment with Z-scores.

## Install

```bash
# LDAK 6+
wget https://raw.githubusercontent.com/dougspeed/LDAK/main/ldak6.3.linux
chmod +x ldak6.3.linux
```

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| LDAK tagging file: build mismatch | Using hg19 tagging on hg38 GWAS sumstats | Tagging files are build-specific; download matched build |
