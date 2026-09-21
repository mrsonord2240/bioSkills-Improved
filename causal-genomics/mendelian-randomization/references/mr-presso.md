# MR-PRESSO

Read when running the global / outlier / distortion tests or extracting outlier SNPs. Verbatim from SKILL.md "MR-PRESSO Outlier Detection".

## MR-PRESSO Outlier Detection

**Goal:** Detect horizontal-pleiotropy outliers, remove them, and test whether the corrected estimate differs from the uncorrected one (distortion test).

**Approach:** Three-step framework: global test (presence of pleiotropy), outlier test (per-SNP), distortion test (effect change after outlier removal).

```r
library(MRPRESSO)

dat_p <- dat[dat$mr_keep, ]  # harmonise_data() keeps dropped palindromes as mr_keep = FALSE rows; MR-PRESSO would fit them
set.seed(42)  # mr_presso()'s global/outlier tests are Monte-Carlo; seed for a reproducible p-value
presso <- mr_presso(
    BetaOutcome = 'beta.outcome', BetaExposure = 'beta.exposure',
    SdOutcome = 'se.outcome', SdExposure = 'se.exposure',
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
    data = dat_p, NbDistribution = 10000,  # >= 10000 for publication-grade p-value precision
    SignifThreshold = 0.05
)

print(presso$`MR-PRESSO results`$`Global Test`)         # any pleiotropy
print(presso$`MR-PRESSO results`$`Distortion Test`)     # change after outlier removal

# Outlier SNPs. `Outlier Test` is NULL unless the global test is significant. Its Pvalue is
# ALREADY Bonferroni-adjusted inside MRPRESSO (raw p x nrow(dat)) and may be a string such as
# "<3e-04", so do not divide the threshold by nrow(dat) again (that double correction flagged 0
# outliers where MRPRESSO's own rule flagged 11, 9 of them planted; checked on MRPRESSO 1.0).
# MRPRESSO's own rule is adjusted P <= SignifThreshold.
ot <- presso$`MR-PRESSO results`$`Outlier Test`
outlier_snps <- character(0)
if (!is.null(ot)) {
    p_adj <- suppressWarnings(as.numeric(sub('^<', '', ot$Pvalue)))
    outlier_snps <- dat_p[rownames(ot)[which(p_adj <= 0.05)], 'SNP']
}
```

`NbDistribution` sets the cost and the floor: the outlier test needs `NbDistribution > nrow(dat) / SignifThreshold` (100 SNPs need > 2000), otherwise MRPRESSO warns "Outlier test unstable" and the outlier P is imprecise. 10000 draws on ~100 SNPs did not finish in 40 minutes on a shared CPU, so run it in the background, and use 3000-5000 while exploring.
