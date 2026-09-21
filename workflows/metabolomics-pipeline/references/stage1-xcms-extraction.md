# Stage 1 -- Feature Extraction (modern xcms 4.x)

**Goal:** Turn centroided mzML into a features-by-samples table, carrying the parameters as part of the result.

**Approach:** Use the `MsExperiment`/`XcmsExperiment` containers with `*Param` objects; align to pooled QC, group AFTER alignment (obiwarp aligns the raw profile directly, so no pre-grouping is needed; the PeakGroups method instead needs group -> align -> regroup because it uses grouped anchor peaks), and treat filled values as imputations. Full parameter rationale (ppm, peakwidth, bw, prefilter) lives in metabolomics/xcms-preprocessing.

```r
library(xcms)
# pd: data.frame, one row per file, with a sample_group column ('QC'/'Control'/'Treatment')
# ionization_mode: 'positive' or 'negative', set from the acquisition method -- carried through
# to Stage 3/5 as defs$mode below (commitment #1, the mode-lock). A mixed-mode study runs this
# whole stage twice, once per mode, and merges the resulting feature tables afterward.
raw <- readMsExperiment(spectraFiles = mzml_files, sampleData = pd)

cwp <- CentWaveParam(ppm = 10, peakwidth = c(2, 20), snthresh = 10,
                     prefilter = c(3, 1000), noise = 1000)   # set from instrument; see xcms-preprocessing
xdata <- findChromPeaks(raw, param = cwp)
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6,
    subset = which(sampleData(xdata)$sample_group == 'QC'), subsetAdjust = 'average'))   # anchor RT alignment on pooled QCs
pdp <- PeakDensityParam(sampleGroups = sampleData(xdata)$sample_group,
                        bw = 5, minFraction = 0.5, binSize = 0.025)
xdata <- groupChromPeaks(xdata, param = pdp)        # group on corrected RT (obiwarp needs no pre-grouping)
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())

feat <- featureValues(xdata, value = 'into')        # features x samples; filled cells are imputations
defs <- featureDefinitions(xdata)                   # mzmed / rtmed per feature, for annotation + mummichog
```
