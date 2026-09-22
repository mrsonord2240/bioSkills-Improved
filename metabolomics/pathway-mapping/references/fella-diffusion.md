<!-- Moved verbatim from SKILL.md (2026-09-21). -->

## Network-Diffusion Enrichment (FELLA)

**Goal:** Return the intermediate enzymes, reactions, and modules that mechanistically link the affected metabolites, not just a ranked pathway list.

**Approach:** Build the KEGG knowledge graph once, then per-analysis map KEGG IDs and run heat diffusion; inspect excluded (unmapped) compounds explicitly.

```r
library(FELLA)

# Build once, reuse. buildGraphFromKEGGREST hits the live KEGG API (slow); cache the DB.
graph <- buildGraphFromKEGGREST(organism = 'hsa')
buildDataFromGraph(keggdata.graph = graph, databaseDir = 'fella_hsa', internalDir = FALSE)
fella.data <- loadKEGGdata(databaseDir = 'fella_hsa', internalDir = FALSE)

cpd_ids <- c('C00022', 'C00186', 'C00158', 'C00042', 'C00122', 'C00041') # KEGG compound IDs only
analysis <- defineCompounds(compounds = cpd_ids, data = fella.data)
getExcluded(analysis)                              # compounds that did not map -- report this

# 'diffusion' is the recommended default; runHypergeom = plain ORA over the graph,
# runPagerank (lowercase r) = directed random walks. The method string is lowercase.
# approx = 'normality' (shown here) is analytic/deterministic, no seed needed. If using
# approx = 'simulation' instead, call set.seed() first -- it resamples niter times and is
# not reproducible run-to-run otherwise.
analysis <- runDiffusion(object = analysis, data = fella.data, approx = 'normality')
results <- generateResultsTable(object = analysis, data = fella.data, method = 'diffusion', threshold = 0.05)
```
