# Duplicate Handling: UMI-aware deduplication

## UMI-Aware Deduplication

For UMI libraries (10x scRNA, ctDNA panels, Twist/IDT/Roche UMI capture), naive markdup destroys information. Use UMI-aware tools:

### umi_tools dedup

Input must be **coordinate-sorted and indexed**.
Pass `--paired` for paired-end libraries: without it the mates are deduplicated independently and the output is silently wrong (5689 vs 2805 records on a paired-end capture BAM). Add `--random-seed=1` for a reproducible output: umi_tools picks among tied reads at random, so without it the count moves by about 1 between runs (5688 or 5689 here).

```bash
# 10x / scRNA -- group by cell barcode + UMI. Check the tags exist first: with absent CB/UB,
# --per-cell writes an EMPTY BAM and still exits 0
samtools view cellranger_possorted.bam | head -1000 | grep -c 'CB:Z:'    # must be > 0
umi_tools dedup --stdin=cellranger_possorted.bam --stdout=dedup.bam \
    --extract-umi-method=tag --umi-tag=UB --cell-tag=CB \
    --per-cell --method=directional --random-seed=1
test "$(samtools view -c dedup.bam)" -gt 0

# Bulk UMI, paired-end (UMI in the RX tag)
samtools sort -o sorted.bam raw.bam && samtools index sorted.bam
umi_tools dedup --stdin=sorted.bam --stdout=dedup.bam --paired \
    --extract-umi-method=tag --umi-tag=RX --method=directional --random-seed=1
```

### fgbio consensus (bulk UMI / ctDNA, best practice for low-VAF detection)

`GroupReadsByUmi` needs the mate mapping-quality (`MQ`) tag on every read (see Common Errors). `samtools fixmate -m` on name-grouped input adds it; alternatively `fgbio SetMateInformation` on queryname-sorted input. Consensus reads are written **unmapped**; re-align them before variant calling. Single-strand and duplex use different grouping strategies and are separate branches:

```bash
# If the UMI is in a separate FASTQ instead of the RX tag, annotate first and use annotated.bam below:
#   fgbio AnnotateBamWithUmis -i raw.bam -f umi.fastq -o annotated.bam
samtools sort -n -o qn.bam raw.bam
samtools fixmate -m qn.bam mated.bam        # or: fgbio SetMateInformation -i qn.bam -o mated.bam

# Single-strand molecular consensus
fgbio GroupReadsByUmi -i mated.bam -o grouped.bam --strategy=adjacency --edits=1 --raw-tag=RX
fgbio CallMolecularConsensusReads -i grouped.bam -o consensus.bam --min-reads=1

# Duplex (xGen-Prism, NEBNext duplex): needs --strategy=paired, which requires RX as two UMIs joined by '-'
# (UMI1-UMI2; a single-UMI RX fails with IllegalArgumentException; single-UMI libraries use the adjacency
# branch above) and writes MI tags with /A /B strand suffixes. CallDuplexConsensusReads on adjacency-grouped reads crashes (StringIndexOutOfBoundsException).
fgbio GroupReadsByUmi -i mated.bam -o grouped_duplex.bam --strategy=paired --edits=1 --raw-tag=RX
fgbio CallDuplexConsensusReads -i grouped_duplex.bam -o duplex.bam --min-reads 1 1 0
```

### Picard UMI-aware marking
```bash
picard UmiAwareMarkDuplicatesWithMateCigar I=coordsort_fixmate.bam O=marked.bam M=metrics.txt \
    UMI_METRICS=umi_metrics.txt UMI_TAG_NAME=RX
```

`--method=directional` is the default and correct -- do not use `--method=unique`, which treats single-base UMI errors as different molecules. `samtools markdup --barcode-tag RX` (UMI/barcode handling added in samtools 1.16) does exact-match UMI grouping; adequate for IDT xGen Duplex but insufficient for single-UMI applications where 1-edit errors are common.
