# Kinase Activity with KSEAapp

## Kinase Activity with KSEAapp

**Goal:** Score kinase activity from the protein-adjusted site fold-changes.

**Approach:** `KSEAapp::KSEA.Scores` merges the site table with a kinase-substrate prior on gene symbol + residue and returns a z-score per kinase with its substrate count `m`. The prior (`KSData`: PhosphoSitePlus + NetworKIN table with `KINASE`, `SUB_GENE`, `SUB_MOD_RSD`, `Source`, `networkin_score` columns) must be the full file from github.com/casecpb/KSEA; the package's `data(KSData)` is an abbreviated demonstration subset. `PX` needs six columns in this exact order, and `FC` is a LINEAR ratio (the function takes log2 itself), treatment over control.

```r
library(KSEAapp)
KSData <- read.csv('PSP&NetworKIN_Kinase_Substrate_Dataset.csv')   # user-supplied prior
pg <- rd('proteinGroups_global.txt')
gene_symbol <- setNames(sub(';.*', '', pg$Gene.names), sub(';.*', '', pg$Protein.IDs))

# log2FC = +/-Inf (site measured in one condition only) gives FC = 0 or Inf, and one such row
# turns every kinase z-score NaN; filter on log2FC, not on the linear FC (0 is finite)
ks <- adjusted[is.finite(adjusted$log2FC), ]
if (nrow(ks) == 0)
  stop('No site has a finite log2FC -- every row of the adjusted table is +/-Inf or NA, i.e. every ',
       'site was quantified in one condition only. KSEA scores fold changes, not presence/absence; ',
       'report those sites as detected-in-one-condition instead.')
PX <- data.frame(
  Protein = sub('_[STY][0-9]+$', '', ks$Protein),
  Gene = gene_symbol[sub('_[STY][0-9]+$', '', ks$Protein)],   # HUGO symbol; merge key with SUB_GENE
  Peptide = rep('NULL', nrow(ks)),   # rep(), not the bare literal: a length-1 value against
                                     # zero-length columns is 'differing number of rows: 0, 1'
  Residue.Both = sub('^.*_', '', ks$Protein),                  # e.g. S473; merge key with SUB_MOD_RSD
  p = ks$adj.pvalue,
  FC = 2^ks$log2FC)   # Treatment/Control because the contrast above is 'Treatment vs Control'
PX <- PX[!is.na(PX$Gene), ]
if (nrow(PX) == 0)
  stop('PX is empty after dropping sites with no gene symbol: proteinGroups_global.txt carried no ',
       '`Gene names` for any tested protein. Check the search FASTA had gene annotation, or map the ',
       'accessions to HUGO symbols yourself before building PX.')

# Check PRIOR COVERAGE before calling. KSEA.Scores merges on SUB_GENE + SUB_MOD_RSD and then
# aggregates, so a prior overlapping this site list in 0 or 1 place dies inside the package with
# `no rows to aggregate` -- a message that says nothing about coverage. Curated priors cover only a
# small fraction of any real site list, so this is the common failure, not an exotic one.
prior <- KSData[grep('PhosphoSitePlus', KSData$Source), ]   # the subset NetworKIN = FALSE will use
covered <- sum(paste(PX$Gene, PX$Residue.Both) %in% paste(prior$SUB_GENE, prior$SUB_MOD_RSD))
cat('sites in PX:', nrow(PX), '| covered by the prior:', covered, '\n')
if (covered < 2)
  stop('The kinase-substrate prior covers ', covered, ' of ', nrow(PX), ' sites, so KSEA has ',
       'nothing to score. Check that KSData is the FULL PhosphoSitePlus+NetworKIN table (not the ',
       'abbreviated data(KSData)), that Gene holds HUGO symbols matching SUB_GENE, and that ',
       'Residue.Both is formatted like SUB_MOD_RSD (S473, not pS473 or Ser473).')

# NetworKIN = FALSE: PhosphoSitePlus-curated pairs only; TRUE adds predictions above NetworKIN.cutoff
# (and then coverage should be counted against that subset instead).
kinase_scores <- KSEA.Scores(KSData, PX, NetworKIN = FALSE, NetworKIN.cutoff = 3)
kinase_scores[order(kinase_scores$z.score), c('Kinase.Gene', 'm', 'z.score', 'FDR')]
```

If the contrast Label reads `Control vs Treatment`, use `FC = 2^(-log2FC)`, or every kinase's sign inverts.
