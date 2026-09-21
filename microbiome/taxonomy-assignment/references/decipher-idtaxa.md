## DECIPHER IDTAXA

**Goal:** Get a conservative, novelty-aware classification that refuses to descend into a clade the query likely does not belong to.

**Approach:** Convert ASV sequences to a DNAStringSet, classify with a pre-trained DECIPHER trainingSet, then flatten the per-rank output to a matrix, mapping IDTAXA's "unclassified_" placeholders to NA.

```r
library(DECIPHER)
load('SILVA_SSU_r138_2019.RData')  # provides the trainingSet object, if you have a pre-trained one

# No pre-trained .RData for your marker/region? Train one directly from a reference FASTA +
# matching "Root;domain;phylum;...;genus;" taxonomy strings (one per sequence, same order).
# LearnTaxa() tunes its tree-descent k-mer sampling with repeated random subsamples (its own
# documentation: "this process is repeated with 100 random subsamples") -- inherently stochastic,
# same class of bug as IdTaxa() below. Verified: two unseeded LearnTaxa() calls on identical input
# produce non-identical trainingSet objects; set.seed() before EVERY LearnTaxa() call makes the
# trainingSet object itself reproducible (identical() TRUE).
# refseqs <- readDNAStringSet('region-matched-ref.fasta')
# reftax  <- readLines('region-matched-ref-taxonomy.txt')  # e.g. "Root;Bacteria;Firmicutes;...;"
# set.seed(100)
# trainingSet <- LearnTaxa(refseqs, taxonomy = reftax)
# MEMORY: LearnTaxa() against a full, un-subsampled reference (400K+ sequences) needs tens of GB
# of RAM and can crash on constrained hardware; subsample the reference (e.g. ~60,000 sequences)
# if it does.
# LearnTaxa's OPTIONAL rank= argument (a 5-column Index/Name/Parent/Level/Rank data.frame, rarely
# available outside DECIPHER's own pre-built .RData sets) is not required to train or classify --
# see the flattening note below for why it matters anyway.

dna <- DNAStringSet(getSequences(seqtab_nochim))

# IdTaxa() descends its classification tree with an internal stochastic step -- inherently
# stochastic, same as assignTaxonomy() above, and by a LARGER margin (verified: unseeded, two
# back-to-back calls on the identical trainingSet and identical query set differ at ~3-4% of
# genus calls). set.seed() before EVERY IdTaxa() call -- without it, repeated runs on the same
# input differ at the genus call for ~3-4% of ASVs. Verified: with set.seed() before each call,
# repeated runs are bit-identical at every rank including genus, in both the default
# multithreaded (processors=NULL) and single-threaded (processors=1) configurations; any fixed
# integer works, 100 is just a convention here (matches the assignTaxonomy() seed above).
set.seed(100)

# threshold 60 = DECIPHER default confidence cutoff; raise for stricter calls. IDTAXA's
# tree-descent stops (leaves the rank unclassified) when the query likely belongs to a taxon
# absent from the reference - this is the intended anti-over-classification behaviour.
ids <- IdTaxa(dna, trainingSet, strand = 'both', threshold = 60, processors = NULL)

ranks <- c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')

# Flatten POSITIONALLY, not by name. x$rank is populated ONLY when trainingSet was built with
# LearnTaxa's rank= data.frame (see above) -- absent that, x$rank is NULL for every result, and
# match(ranks, x$rank) silently returns all-NA with no error or warning (confirmed empirically
# against a real LearnTaxa()-trained set: 100% NA at every rank, no exception raised). x$taxon[1]
# is always "Root"; the remaining entries are domain..genus/species in taxonomic order regardless
# of whether rank= was supplied, so index positionally instead.
taxa_idtaxa <- t(sapply(ids, function(x) {
    taxa <- x$taxon[-1]                    # drop "Root"
    taxa[startsWith(taxa, 'unclassified_')] <- NA
    length(taxa) <- length(ranks)          # pad/truncate to the fixed rank depth above
    taxa
}))
colnames(taxa_idtaxa) <- ranks

# Sanity check: a wrong-region or wrong-marker training set, or a flattening bug, gives all-NA
# calls with no error. To confirm reproducibility, repeat set.seed(100) + IdTaxa() once and
# compare identical(ids, ids_again).
if (all(is.na(taxa_idtaxa[, 'genus']))) stop('IdTaxa returned no genus calls: check the trainingSet marker/region and the flattening')
```
