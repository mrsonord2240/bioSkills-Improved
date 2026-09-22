# Tool installation

Read when a package is missing or `install.packages` fails. Verbatim from SKILL.md "Tool Installation Notes".

## Tool Installation Notes

```r
# CRAN-stable
install.packages(c('remotes', 'MendelianRandomization', 'MVMR', 'coloc', 'simex'))

# GitHub-only or recently archived
remotes::install_github('MRCIEU/TwoSampleMR')          # primary orchestrator
remotes::install_github('MRCIEU/ieugwasr')             # OpenGWAS client + local clumping
remotes::install_github('rondolab/MR-PRESSO')          # never on CRAN
remotes::install_github('qingyuanzhao/mr.raps')        # CRAN-archived 2025-03-01
remotes::install_github('jean997/cause')               # depends on mixsqp; suggests Rfast
remotes::install_github('cnfoley/mrclust')             # heterogeneity clusters
remotes::install_github('LizaDarrous/lhcMR')           # bidirectional + heritable confounder
remotes::install_github('HDTian/DRMR')                 # doubly-ranked stratification
remotes::install_github('n-mounier/MRlap')             # joint overlap + winner's-curse + weak-IV correction
```

`TwoSampleMR::mr_raps()` is a thin wrapper that calls `mr.raps::mr.raps()` under the hood; the GitHub `mr.raps` install above is therefore required. The `MendelianRandomization` package does NOT export `mr_raps()` (verify with `ls('package:MendelianRandomization')`); only TwoSampleMR offers a MR-RAPS entry point. For local clumping, install plink2 (https://www.cog-genomics.org/plink/2.0/) and download a 1KG EUR (or matched-ancestry) reference bfile (prebuilt at https://mrcieu.github.io/ieugwasr/).
