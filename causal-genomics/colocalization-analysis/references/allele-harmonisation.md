# Allele harmonisation

## Allele Harmonisation (Critical Pre-Step)

Mismatched effect alleles silently invert signs of betas, collapsing PP.H4 into PP.H3. Required steps before coloc:

1. Merge GWAS and eQTL summary stats by SNP ID (rsID or chr:pos:ref:alt).
2. Mark SNP-pairs as `same` (A1/A2 match) or `flip` (A1/A2 swap); for non-palindromic pairs reported on opposite strands (e.g. A/G vs T/C), complement dataset 2's alleles first, then classify; drop SNPs that match neither.
3. For `flip` rows, negate the second dataset's beta (and swap A1/A2).
4. Drop palindromic SNPs (A/T or C/G) at MAF > 0.42; their strand cannot be inferred from coding alone (TwoSampleMR `harmonise_data` standard cutoff).
5. Verify genome build alignment (hg19 vs hg38 must match; lift over if not).

```r
harmonise <- function(df1, df2) {
    comp <- function(a) c(A='T', T='A', C='G', G='C')[a]
    m <- merge(df1, df2, by='SNP', suffixes=c('.1','.2'))
    palindromic <- (m$A1.1 %in% c('A','T') & m$A2.1 %in% c('A','T')) |
                   (m$A1.1 %in% c('C','G') & m$A2.1 %in% c('C','G'))
    # Non-palindromic strand mismatch: complement dataset 2's alleles, then treat as same/flip
    strand <- !palindromic & m$A1.1 == comp(m$A1.2) & m$A2.1 == comp(m$A2.2) |
              !palindromic & m$A1.1 == comp(m$A2.2) & m$A2.1 == comp(m$A1.2)
    strand[is.na(strand)] <- FALSE
    a1 <- ifelse(strand, comp(m$A1.2), m$A1.2); a2 <- ifelse(strand, comp(m$A2.2), m$A2.2)
    m$A1.2 <- unname(a1); m$A2.2 <- unname(a2)
    same <- m$A1.1 == m$A1.2 & m$A2.1 == m$A2.2
    flip <- m$A1.1 == m$A2.2 & m$A2.1 == m$A1.2
    m$BETA.2[flip] <- -m$BETA.2[flip]
    m$MAF.2[flip] <- 1 - m$MAF.2[flip]
    keep <- (same | flip) & !(palindromic & m$MAF.1 > 0.42)
    m[keep, ]
}
```

Harmonisation pitfalls to watch for:

- **Allele coding mismatch.** GWAS may report effect allele as A1 while eQTL reports it as A2. Always check both and flip betas where needed.
- **Build mismatch.** hg19 GWAS coords + hg38 eQTL coords silently merge on rsID but break on chr:pos. Lift over with `rtracklayer::liftOver` or CrossMap before merging.
- **Strand mismatch.** Non-palindromic SNPs coded on opposite strands (A/G vs T/C) are resolved by the complement step in `harmonise()`. Palindromic SNPs cannot be resolved this way (next bullet).
- **Palindromic SNPs at high MAF.** A/T and C/G SNPs at MAF > 0.42 cannot be unambiguously strand-resolved; drop them or resolve with reference-panel MAF.
- **Multi-allelic SNPs.** Many summary stats collapse multi-allelic loci by keeping only the most-frequent alt; if datasets pick different alts, harmonisation drops the SNP. Split on chr:pos:ref:alt as a unique key.
- **rsID dependence.** rsID can be remapped across dbSNP builds (e.g. merge of two rsIDs into one). Prefer chr:pos:ref:alt keys for cross-study merges.
