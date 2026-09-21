# FINEMAP CLI

## FINEMAP CLI Pattern

**Goal:** Independent confirmation via shotgun stochastic search.

**Approach:** Build .z, .ld, and master files; run FINEMAP with `--sss` and parse the .snp and .cred outputs.

```bash
# .z file format: snp chromosome position allele1 allele2 maf beta se
# .ld file: square LD matrix, space-separated, no header

cat > locus.master <<'EOF'
z;ld;snp;config;cred;log;n_samples
locus.z;locus.ld;locus.snp;locus.config;locus.cred;locus.log;500000
EOF

finemap --sss \
    --in-files locus.master \
    --n-causal-snps 5 \
    --prob-tol 0.001 \
    --n-iterations 100000 \
    --n-convergence 5000

# Parse:
# locus.snp -> per-variant prob (PIP), log10bf
# locus.cred -> credible sets at increasing causal counts
# locus.config -> top configurations
```

FINEMAP and SuSiE agree when sparsity holds; disagreement often reveals non-sparse loci that need SuSiE-inf.
