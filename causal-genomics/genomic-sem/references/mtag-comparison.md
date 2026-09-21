# MTAG comparison and reconciliation (reference for genomic-sem SKILL.md)

Running MTAG beside GenomicSEM, and what to do when the two disagree. Read it when `SKILL.md` "Reference Files" points here. The section names below refer to `SKILL.md`.

## MTAG Comparison

**Goal:** Cross-check GenomicSEM common-factor results against MTAG per-trait shrunk z-scores.

**Approach:** Run MTAG CLI on the same input sumstats; compare top hits with GenomicSEM factor hits. Report MaxFDR.

```bash
# MTAG CLI (Python)
python mtag.py \
    --sumstats trait1.txt,trait2.txt,trait3.txt \
    --n_min 0 \
    --out mtag_results
# MTAG uses the signed Z by default; --use_beta_se was disabled upstream (raises a
# RuntimeError since Dec 2021 due to beta-se bugs), so supply a Z column and omit it.

# Check MaxFDR per trait
grep -iE 'max ?fdr' mtag_results.log   # matches both the section header and the 'Max FDR of Trait' value lines
# Each per-trait MTAG file: mtag_results_trait_<k>.txt
```

If MaxFDR > 0.05 for any trait, MTAG results for that trait are unreliable; GenomicSEM with Q_SNP filtering is the more defensible report.

## Reconciliation: When GenomicSEM and MTAG Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| GenomicSEM factor SNP sig, MTAG sig for all traits | Genuine common-factor SNP | Report; high confidence |
| GenomicSEM factor SNP sig, MTAG sig in only 1 trait | Q_SNP heterogeneity likely; one-trait-dominant | Check Q_SNP; if sig, this is NOT a factor SNP |
| MTAG sig, GenomicSEM factor null, Q_SNP sig | Trait-specific SNP captured by MTAG shrinkage | Report as trait-specific, not common-factor |
| Both null but per-trait univariate sig | Power loss from multivariate parameterization | Re-check sample overlap V matrix |
| GenomicSEM and MTAG both sig but opposite direction | Sample-overlap mis-specification OR sign error in munging | Re-munge with same allele convention; re-run `ldsc()` |
| MTAG MaxFDR > 5%, GenomicSEM with Q_SNP works | MTAG assumption violated | Prefer GenomicSEM as primary |
| One-trait GWAS sig but common-factor not | Trait-specific architecture | Don't force into common-factor frame |

**Operational rule for publication:** A common-factor SNP claim requires (1) factor p < 5e-8, (2) Q_SNP p > 0.05 / N_factor_SNPs (non-heterogeneous), and (3) replication in an independent set of traits or cohorts. Trait-specific SNPs from MTAG require MaxFDR < 5% for the trait. Reporting only the factor effect without Q_SNP is the most common reviewer-flagged error.
