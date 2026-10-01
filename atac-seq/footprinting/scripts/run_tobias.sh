#!/bin/bash
# Tested with TOBIAS 0.17.5, samtools 1.19, bedtools 2.31 | Verify flags against `TOBIAS <tool> --help` if the version differs
# TOBIAS three-step ATAC-seq footprinting: bias correction, footprint scoring, differential bound/unbound calls,
# a CTCF positive-control plot per condition, and a ranked differential summary.
#
# Usage: run_tobias.sh cond1.bam cond2.bam peaks.bed genome.fa blacklist.bed motifs.pfm [outdir] [cores]
# Environment: ALLOW_NO_CTCF=1  continue (exit 0) when the motif file yields no CTCF motif; default is to exit 3
#              PMAX=0.05        p-value ceiling for the differential summary

set -euo pipefail

COND1_BAM=${1:-cond1.dedup.nochrM.bam}
COND2_BAM=${2:-cond2.dedup.nochrM.bam}
PEAKS=${3:-consensus_peaks.bed}                  # Must be the same peakset for both conditions
GENOME=${4:-hg38.fa}
BLACKLIST=${5:-hg38-blacklist.v2.bed}
MOTIFS=${6:-JASPAR2024_CORE_vertebrates.pfm}
OUTDIR=${7:-tobias_out}
CORES=${8:-16}
PMAX=${PMAX:-0.05}

for f in "$COND1_BAM" "$COND2_BAM" "$PEAKS" "$GENOME" "$BLACKLIST" "$MOTIFS"; do
    [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 2; }
done

# Print the single file matching a glob pattern inside a directory; fail on zero or several matches.
one_match() {
    local m=( "$1"/$2 )
    if [ "${#m[@]}" -ne 1 ] || [ ! -f "${m[0]}" ]; then
        echo "ERROR: expected exactly one file matching $2 in $1, found ${#m[@]}" >&2
        return 1
    fi
    printf '%s
' "${m[0]}"
}

mkdir -p "$OUTDIR"/{cond1,cond2,bindetect,validation}

# Step 1: Bias correction per condition (matched peakset and blacklist)
for cond in cond1 cond2; do
    [ "$cond" = cond1 ] && BAM=$COND1_BAM || BAM=$COND2_BAM
    echo "=== ATACorrect: $cond ==="
    TOBIAS ATACorrect \
        --bam "$BAM" --genome "$GENOME" \
        --peaks "$PEAKS" --blacklist "$BLACKLIST" \
        --outdir "$OUTDIR/$cond" \
        --cores "$CORES"
    one_match "$OUTDIR/$cond" "*_corrected.bw" > /dev/null
done

# Step 2: Per-base footprint scores
for cond in cond1 cond2; do
    echo "=== ScoreBigwig: $cond ==="
    CORRECTED=$(one_match "$OUTDIR/$cond" "*_corrected.bw")
    TOBIAS ScoreBigwig \
        --signal "$CORRECTED" \
        --regions "$PEAKS" \
        --output "$OUTDIR/$cond/${cond}_footprints.bw" \
        --cores "$CORES"
    [ -s "$OUTDIR/$cond/${cond}_footprints.bw" ] || { echo "ERROR: ScoreBigwig wrote no $cond footprint track" >&2; exit 1; }
done

# Step 3: Differential bound/unbound classification across conditions
echo "=== BINDetect (differential) ==="
TOBIAS BINDetect \
    --motifs "$MOTIFS" \
    --signals "$OUTDIR/cond1/cond1_footprints.bw" "$OUTDIR/cond2/cond2_footprints.bw" \
    --genome "$GENOME" --peaks "$PEAKS" \
    --outdir "$OUTDIR/bindetect" \
    --cond-names cond1 cond2 \
    --cores "$CORES"
RESULTS="$OUTDIR/bindetect/bindetect_results.txt"
[ -s "$RESULTS" ] || { echo "ERROR: BINDetect wrote no $RESULTS" >&2; exit 1; }

# Positive control: CTCF aggregate at all, bound and unbound sites, uncorrected vs corrected signal, per condition.
# Each condition is plotted at its own bound sites. Motif directory names follow the JASPAR name (CTCF_<ID>).
CTCF_ALL=( "$OUTDIR"/bindetect/CTCF_*/beds/CTCF_*_all.bed )
NO_CTCF=0
if [ ! -f "${CTCF_ALL[0]}" ]; then
    NO_CTCF=1
    echo "WARNING: no CTCF motif in $MOTIFS produced beds; the positive-control plot was NOT made." >&2
    echo "         Add a CTCF PFM or supply another positive control before trusting the calls." >&2
else
    CTCF_DIR=$(dirname "$(dirname "${CTCF_ALL[0]}")")
    CTCF_ID=$(basename "$CTCF_DIR")               # first CTCF variant in sort order (MA0139.2 in JASPAR 2024)
    for cond in cond1 cond2; do
        TOBIAS PlotAggregate \
            --TFBS "$CTCF_DIR/beds/${CTCF_ID}_all.bed" \
                   "$CTCF_DIR/beds/${CTCF_ID}_${cond}_bound.bed" \
                   "$CTCF_DIR/beds/${CTCF_ID}_${cond}_unbound.bed" \
            --TFBS-labels all bound unbound \
            --signals "$(one_match "$OUTDIR/$cond" "*_uncorrected.bw")" "$(one_match "$OUTDIR/$cond" "*_corrected.bw")" \
            --signal-labels uncorrected corrected \
            --output "$OUTDIR/validation/ctcf_${cond}_aggregate.pdf" \
            --output-txt "$OUTDIR/validation/ctcf_${cond}_aggregate.txt" \
            --share-y both --plot-boundaries
        echo "Validation plot: $OUTDIR/validation/ctcf_${cond}_aggregate.pdf ($CTCF_ID)"
    done
    echo "Expect a central dip with flanking peaks at bound sites and a flatter profile at unbound sites."
    echo "A dip at bound sites alone does not prove the bias correction worked; compare the uncorrected panel."
fi

# Top differential motifs: |change| descending among p <= PMAX. Motif variants (for example three CTCF IDs) are separate rows.
echo "=== Top differential motifs (p <= $PMAX, ranked by |cond1_cond2_change|) ==="
awk -F'\t' -v pmax="$PMAX" 'BEGIN {OFS="\t"}
    NR==1 {for(i=1;i<=NF;i++) col[$i]=i
           if(!("output_prefix" in col && "cond1_cond2_change" in col && "cond1_cond2_pvalue" in col)) {print "ERROR: unexpected BINDetect columns" > "/dev/stderr"; exit 1}
           next}
    {ch=$(col["cond1_cond2_change"]); p=$(col["cond1_cond2_pvalue"])
     if (p+0 <= pmax) print $(col["output_prefix"]), ch, p, (ch<0 ? -ch : ch)}' "$RESULTS" | \
    sort -t$'\t' -k4,4gr | cut -f1-3 | head -20

if [ "$NO_CTCF" = 1 ] && [ "${ALLOW_NO_CTCF:-0}" != 1 ]; then
    echo "ERROR: no CTCF positive control (set ALLOW_NO_CTCF=1 to accept this)." >&2
    exit 3
fi
