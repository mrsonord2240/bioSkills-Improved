# Alternative Front End -- MS-DIAL

When peak detection happens in the MS-DIAL GUI/console (MS2Dec deconvolution, GC-EI, DIA/SWATH), import the alignment-result table and enter the pipeline at Stage 2. The framing is unchanged: the imported table is still a parameterized hypothesis. See metabolomics/msdial-preprocessing for the export-parsing details, then continue with normalization-qc onward.

The block below turns the export into the objects Stage 2 expects (`feat`, `defs`, `sample_class`, `injection_order`, `batch_id`); run Stage 2 from `fm <- feat` on. Checked on a real MSDIALCUI 5.5.260820 export, whose 4 header rows (`Class`, `File type`, `Injection order`, `Batch ID`) precede the column header and whose `Class` cell marks where the per-sample columns start.

```r
export_file <- Sys.glob('AlignResult-*.mdalign')[1]
hdr <- strsplit(readLines(export_file, n = 4), '\t', fixed = TRUE)   # class / file type / order / batch
msdial <- read.csv(export_file, sep = '\t', skip = 4, check.names = FALSE)
s_idx <- (which(hdr[[1]] == 'Class') + 1):length(hdr[[1]])          # sample columns follow the "Class" cell
file_type <- hdr[[2]][s_idx]
use <- file_type %in% c('Sample', 'QC')      # Blank/Standard injections leave the matrix (blank filter: normalization-qc)

feat <- as.matrix(msdial[, colnames(msdial)[s_idx][use]]); storage.mode(feat) <- 'numeric'
feat[feat == 0] <- NA                        # MS-DIAL writes not-detected as 0; Stage 2's detection filter counts NA
rownames(feat) <- paste0('FT', msdial[['Alignment ID']])
defs <- data.frame(mzmed = msdial[['Average Mz']], rtmed = msdial[['Average Rt(min)']] * 60,  # seconds, as xcms
                   row.names = rownames(feat))
sample_class <- ifelse(file_type[use] == 'QC', 'QC', hdr[[1]][s_idx][use])   # pooled QCs must be labelled 'QC'
injection_order <- as.integer(hdr[[3]][s_idx][use])
batch_id <- as.integer(hdr[[4]][s_idx][use])
stopifnot(ncol(feat) == length(sample_class), !anyNA(injection_order), !anyNA(batch_id))
```
