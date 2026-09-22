<!-- Moved verbatim from SKILL.md (2026-09-21). -->

## Mummichog / PSEA on a Raw m/z Peak Table

**Goal:** Predict perturbed pathway activity from an untargeted LC-MS feature table when no compound identities exist.

**Approach:** Declare instrument ppm and ionization mode, load the FULL feature table (m/z + p-value + t-score, optionally RT), set the query-defining p-cutoff, and run PSEA whose permutation null is sampled from R_all.

```r
current.msg <- character(0); err.vec <- character(0)  # required -- see Version Compatibility
library(MetaboAnalystR)
set.seed(123)  # PerformPSEA's permNum resampling is not reproducible run-to-run otherwise

mSet <- InitDataObjects('mass_all', 'mummichog', FALSE)
mSet <- SetPeakFormat(mSet, 'mpt')               # 'mpt' = m/z, p-value, t-score; 'mprt' adds RT (use with 'v2')

# ppm and ionization mode are chemistry-specific and mandatory; pos and neg use
# entirely different adduct tables. Mixed data needs a per-feature mode column.
mSet <- UpdateInstrumentParameters(mSet, 5.0, 'negative')

# CRITICAL: peaks.txt must be the ENTIRE feature table, not just significant peaks.
# The permutation null draws random feature lists from this file (R_all); supplying
# only significant features pre-enriches the pool and makes everything significant.
mSet <- Read.PeakListData(mSet, 'peaks.txt')
mSet <- SanityCheckMummichogData(mSet)          # also merges duplicate m/z-matched features automatically

mSet <- SetPeakEnrichMethod(mSet, 'mum', 'v2')   # 'mum'|'gsea'|'integ'; 'v2' uses RT/empirical compounds
mSet <- SetMummichogPval(mSet, 0.2)              # query-defining cutoff; default is NOT 0.05 -- document it
mSet <- PerformPSEA(mSet, 'hsa_mfn', 'current', permNum = 1000) # library string encodes organism+network

psea <- mSet$mummi.resmat                         # predicted-active pathways; NOT a metabolite ID list
```

`Read.PeakListData`/`SanityCheckMummichogData` log lines such as `A total of 11 of duplicates were
merged` (four passes on the audit's 1500-feature table): duplicate m/z-matched features are merged
before enrichment, so the feature count PSEA reports can be lower than the input row count.
