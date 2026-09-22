# SPIA and graphite: signed-topology perturbation (KEGG third generation)

Read when scoring signed pathway perturbation on SIGNALING maps. Inputs `de`, `universe` and `org.Hs.eg.db` come from "Prepare the Gene IDs" in SKILL.md.

## Run Signed-Topology Perturbation (SPIA) -- the Third Generation

**Goal:** Score how perturbed each SIGNALING pathway is given both the over-representation of DE genes and the propagation of their fold-changes through the signed wiring.

**Approach:** SPIA combines pNDE (the classical over-representation evidence) with pPERT (the probability of the observed total accumulated perturbation tA, computed by propagating log2 fold-changes through KGML activation/inhibition edges) into a single global pG, then FDR-corrects it. It needs a NAMED vector of DE fold-changes plus the universe, and is defined only for signaling maps. Two routes: `spia()` reads SPIA's bundled `hsaSPIA` KEGG topology (the default for direction calls); graphite harmonizes node IDs, resolves complexes/families, removes compounds and reads current KEGG (or Reactome) topology, but its perturbation direction can differ (see the caveat below the code).

```r
library(SPIA)
sig <- de[de$padj < 0.05, ]   # DE genes only
map <- bitr(sig$gene, 'SYMBOL', 'ENTREZID', org.Hs.eg.db)   # bitr drops/many-to-one: MERGE, never assign as names
de_vec <- setNames(sig$log2FoldChange[match(map$SYMBOL, sig$gene)], map$ENTREZID)
de_vec <- de_vec[!duplicated(names(de_vec))]
set.seed(123)   # SPIA's pPERT is a stochastic bootstrap; fix the seed before spia()
res <- spia(de=de_vec, all=universe, organism='hsa', nB=2000, plots=FALSE)   # nB=2000 bootstraps for pPERT
# output cols: Name, ID, pSize, NDE, pNDE, tA, pPERT, pG, pGFdr, pGFWER, Status, KEGGLINK
# Status reports inferred Activated / Inhibited from the sign of tA

# graphite route (current KEGG topology instead of SPIA's bundled snapshot; works on Reactome too)
library(graphite)
db <- pathways('hsapiens', 'kegg')
db <- convertIdentifiers(db, 'ENTREZID')
# runSPIA checks `datasetName(pathwaySetName) %in% dir()`, and bare dir() lists only the
# CURRENT WORKING DIRECTORY's filenames -- an absolute/tempdir() pathwaySetName can never
# match, so prepareSPIA/runSPIA must both run with a RELATIVE name from a matching setwd().
# convertIdentifiers() also prefixes graphite's node IDs ('ENTREZID:1017'), so de_vec/all
# need the same prefix or every ID join returns 0 rows even once the path bug is worked
# around. Confirmed against installed graphite 1.52.0 and current Bioconductor-release
# graphite 1.56.0 source.
de_vec_gr  <- setNames(de_vec, paste0('ENTREZID:', names(de_vec)))
universe_gr <- paste0('ENTREZID:', universe)
owd <- getwd(); setwd(tempdir())
prepareSPIA(db, 'kegg_hsa_spia')              # writes kegg_hsa_spiaSPIA.RData into tempdir()
set.seed(123)   # graphite's runSPIA bootstraps pPERT the same way spia() does
gr <- runSPIA(de=de_vec_gr, all=universe_gr, 'kegg_hsa_spia')
setwd(owd)
```

**The two routes are complementary evidence, not interchangeable.** Same `de_vec`, universe and seed, they score different topologies (bundled `hsaSPIA`: 139 pathways, an older snapshot; graphite: 319 graphs from the live KEGG conversion, e.g. 233 binding/association and 269 inhibition edges in Cell cycle) and disagree on direction for a meaningful fraction of pathways. On the audit's synthetic data (nB=50; 99 pathways scored by both) tA correlated only r=0.60, 14 pathways had opposite-sign tA, 18 Activated/Inhibited calls differed (further pathways had tA=0 in one route), and Cell cycle (planted UP) was Activated in `spia()` but Inhibited in graphite (re-run at nB=100: same split; SPIA 2.58.0, graphite 1.52.0). Report Activated/Inhibited only where both routes agree, state which topology produced a single-route call, and prefer `spia()` when only one is run; use graphite when current KEGG or Reactome topology is required.

SPIA aborts if more than ~1% of the DE IDs are absent from `all`, so build the universe from the same ID space. The standalone SPIA package also ships a frozen `hsaSPIA` data object that is an OLDER snapshot than a live enrichKEGG query - do not mix the two in one comparison.
