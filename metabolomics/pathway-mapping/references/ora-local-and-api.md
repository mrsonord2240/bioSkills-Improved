<!-- Moved verbatim from SKILL.md "ORA on an Identified Compound List" (2026-09-21). Run the mapping block in SKILL.md first: it produces `kegg_ids`. -->

### Local-Only ORA (default: no remote call on user data, verified background correction)

Computes the hypergeometric ORA locally from KEGGREST's public pathway-to-compound table (generic
reference data, not user data). Checked on KEGGREST 1.46.0: ~3 s to fetch and build the table. On the
audit's synthetic 12-compound TCA-cycle input the Citrate cycle (hsa00020) p-value moved from
6.1e-19 (all-of-KEGG background, n=6709) to 9.4e-11 (320-ID assay-coverage background; re-run 2026-09-21) --
less significant with the correct, smaller background, as the theory predicts.

```r
library(KEGGREST)

local_kegg_ora <- function(hit_kegg_ids, universe_kegg_ids, min_hits = 2) {
  links <- keggLink('pathway', 'compound')             # public reference table; no user data sent
  cpd_ids  <- sub('^cpd:', '', names(links))
  path_ids <- sub('^path:map', 'hsa', unname(links))    # generic map#### -> organism-specific hsa####
  pw2cpd <- lapply(split(cpd_ids, path_ids), unique)

  universe_kegg_ids <- unique(universe_kegg_ids)
  pw2cpd <- lapply(pw2cpd, function(x) intersect(x, universe_kegg_ids))
  pw2cpd <- pw2cpd[lengths(pw2cpd) > 0]
  hits <- intersect(hit_kegg_ids, universe_kegg_ids)
  N <- length(universe_kegg_ids); k <- length(hits)

  out <- data.frame(
    pathway = names(pw2cpd),
    total   = lengths(pw2cpd),
    hits    = vapply(pw2cpd, function(s) length(intersect(s, hits)), integer(1))
  )
  out <- out[out$hits >= min_hits, ]
  out$p.value <- mapply(function(m, h) phyper(h - 1, m, N - m, k, lower.tail = FALSE), out$total, out$hits)
  out$fdr <- p.adjust(out$p.value, method = 'BH')
  out[order(out$p.value), ]
}

# reference_ids: the assay-coverage background -- one KEGG compound ID per line
reference_ids <- readLines('reference_metabolome.txt')
reference_ids <- reference_ids[nzchar(reference_ids)]
ora_local <- local_kegg_ora(kegg_ids, reference_ids)   # kegg_ids from the mapping block above
```

### MetaboAnalystR API path (alternative; sends the compound list off-machine)

Use only if the remote call disclosed in Version Compatibility is acceptable. Run the mapping block above first.

```r
mSet <- SetKEGG.PathLib(mSet, 'hsa', 'current')

# SetMetabolomeFilter(mSet, TRUE) alone does NOT restrict the background --
# Setup.KEGGReferenceMetabolome() must be called first to load the reference file into
# mSet$dataSet$metabo.filter.kegg, or the filter silently has no effect. FALSE uses the
# whole library (all of KEGG) -- the inflated default that manufactures false positives.
mSet <- SetMetabolomeFilter(mSet, TRUE)
mSet <- Setup.KEGGReferenceMetabolome(mSet, 'reference_metabolome.txt')  # one KEGG ID per line

mSet <- CalculateOraScore(mSet, 'rbc', 'hyperg') # node-importance 'rbc'|'dgr'; test 'hyperg'|'fisher'
# Checked on MetaboAnalystR 4.3.0: the server has been observed to reject the FILTERED
# request outright (CalculateOraScore returns 0; current.msg == "Failed to connect to
# the API Server!"), even with a correctly-matched reference file, while the unfiltered
# (FALSE) call to the same endpoint succeeds. If this happens, use Local-Only ORA above.
if (is.numeric(mSet)) {
  cat('ORA failed:', paste(current.msg, collapse = ' | '), '\n')
} else {
  ora <- as.data.frame(mSet$analSet$ora.mat)   # columns include Raw p, FDR, Impact, Hits, Total
}
```
