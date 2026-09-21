# FINEMAP CLI

## FINEMAP CLI Pattern

**Goal:** Independent confirmation via shotgun stochastic search.

**Approach:** Build .z, .ld, and master files; run FINEMAP with `--sss` and parse the .snp and .cred outputs.

Full pipeline (locus extraction, `plink --r square` LD, `.z` and master-file construction, run, report): `examples/finemap_pipeline.sh` (`.z` columns: rsid chromosome position allele1 allele2 maf beta se; `.ld` is a square, space-separated, headerless signed-r matrix; master header `z;ld;snp;config;cred;log;n_samples`). Tighter search than the example's defaults:

```bash
finemap --sss --in-files locus.master --n-causal-snps 5 --prob-tol 0.001 --n-iterations 100000 --n-convergence 5000
```

Outputs: `locus.snp` (per-variant PIP `prob`, `log10bf`), `locus.cred` (credible sets at increasing causal counts), `locus.config` (top configurations).

FINEMAP and SuSiE agree when sparsity holds; disagreement often reveals non-sparse loci that need SuSiE-inf.
