# SuSiEx Cross-Ancestry Fine-Mapping

## Cross-Ancestry Fine-Mapping with SuSiEx

**Goal:** Jointly fine-map a locus across multiple ancestries assuming shared causal variants but population-specific LD.

**Approach:** Per-ancestry summary statistics + per-ancestry LD reference; SuSiEx runs a joint SuSiE model with population-specific R matrices. SuSiEx assigns populations by the ORDER of the comma-separated `--sst_file`/`--n_gwas`/`--ref_file`/`--ld_file` lists (there is no `--pop` flag); keep all four lists in the same population order.

Full command with every required column flag (`--chr_col` ... `--pval_col`, one value per population), `--ld_file`, `--level` and output layout: `examples/susiex_multiancestry.sh` (also the one-line form in `SKILL.md`).

The output includes per-population PIPs and a joint credible set. Credible sets from SuSiEx are typically 2-5x smaller than EUR-only susie_rss when AFR is included, because AFR shorter LD blocks resolve EUR-tagged regions.
