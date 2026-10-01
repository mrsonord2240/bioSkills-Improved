# Footprinting prerequisites and request examples

The core Skill describes the workflow; [`method-reference.md`](method-reference.md) holds method details and caveats.

## Environments

The tools pin incompatible Pythons (rgt and pydnase need Python 3.7-era stacks, scPrinter needs 3.11 and pins), so use one environment per tool. Installing them together resolves TOBIAS 0.13.3 on Python 3.7, below the tested version, and `bioconda` without `conda-forge` cannot solve `tobias` (missing `adjusttext`). These are the commands the tested environments were built from (micromamba shown; conda or mamba work alike):

```bash
# TOBIAS 0.17.5 (also samtools, bedtools, deepTools alignmentSieve)
micromamba create -n footprint -c conda-forge -c bioconda python=3.10 pip pybigwig pysam samtools=1.19 bedtools numpy deeptools
micromamba run -n footprint pip install tobias

# HINT-ATAC (rgt 1.0.2) and Wellington (pyDNase 0.3.0), each alone
micromamba create -n footprint-rgt -c conda-forge -c bioconda rgt
micromamba create -n footprint-pydnase -c conda-forge -c bioconda pydnase

# scPrinter 1.2.0 (GitHub tag, not on PyPI): torch build for your CUDA, then pins
micromamba create -n footprint-scprinter -c conda-forge -c bioconda python=3.11 pip git bedtools
pip install torch --index-url https://download.pytorch.org/whl/cu128     # tested: torch 2.11.0+cu128
pip install "git+https://github.com/buenrostrolab/scPrinter@v1.2.0" "tangermeme==0.4.4" "snapatac2==2.8.0" ema_pytorch
```

The scPrinter pins exist because v1.2.0 imports `tangermeme.tools.tomtom` (absent in tangermeme 1.5.0) and `import_fragments` needs `snapatac2.pp.import_data` (absent in 2.9.0); snapatac2 2.8.0 downgrades numpy, pandas, and anndata and pip warns about zarr and shap, but the tested path ran. The pinned command was re-run in a fresh environment (about 30 min); it builds macs3, MOODS-python, and sorted_nearest from source, so it needs a C/C++ compiler and git. Set `SCPRINTER_DATA` to a writable directory; the first use downloads pretrained models and needs internet.

### RGT data for HINT-ATAC

The bioconda `rgt` package ships no data directory, and `rgt-hint` fails without `data.config`. Reinstall RGT from PyPI source with `RGTDATA` set once (the bioconda copy makes a plain `pip install` a no-op, so `--force-reinstall --no-cache-dir` is required; `setup.py` then writes `data.config` and copies the HMMs and bias tables), then point the hg38 entry at your own genome and annotation (or run `python setupGenomicData.py --hg38` to download them):

```bash
export RGTDATA=$HOME/rgtdata
pip install --no-deps --force-reinstall --no-cache-dir rgt==1.0.2   # in the footprint-rgt environment; writes $RGTDATA/data.config
cd "$RGTDATA" && python setupGenomicData.py --hg38 --hg38-genome-path hg38.fa --hg38-gtf-path gencode.gtf
```

### Fragment file for scPrinter

Raw, unshifted fragments of properly paired reads, coordinate-sorted, bgzipped and indexed with the `footprint` environment's samtools and htslib (scPrinter applies the +4/-5 shift itself):

```bash
samtools view -f 2 sample.bam | awk 'BEGIN{OFS="\t"} $9>0 && $9<1000 {print $3,$4-1,$4-1+$9,"sample"}' \
    | sort -k1,1 -k2,2n -k3,3n > frags.tsv && bgzip frags.tsv && tabix -p bed frags.tsv.gz
```

## Inputs

- Deduplicated, MAPQ-filtered, chrM-stripped BAM (at least 50M nuclear reads recommended).
- Consensus peakset (BED or narrowPeak).
- Reference genome FASTA matching the BAM build.
- Assembly-matched blacklist.
- Motif database: JASPAR 2024 CORE vertebrates PFM (default) or HOCOMOCO v12.

```bash
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
mv JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt JASPAR2024_CORE_vertebrates.pfm   # default motif name in scripts/run_tobias.sh
```

## Request examples

- "Run the TOBIAS three-step pipeline (ATACorrect, ScoreBigwig, BINDetect) on two conditions, with `--cond-names treated control`."
- "Plot the aggregate footprint at JASPAR CTCF MA0139.2 sites with TOBIAS PlotAggregate and confirm a clean dip at bound sites."
- "Run scPrinter for multi-scale footprints on this bulk library and compare CTCF bound and unbound sites."
- "Run TOBIAS and HINT-ATAC; report the bound-versus-unbound site overlap of the HINT footprints as the concordance measure."
- "Filter to fragments under 100 bp, run TOBIAS on the NFR BAM, and compare footprint sharpness with the full-fragment BAM."
- "Cross-validate predicted FOXA1 footprints against published ChIP-seq peaks."

Single-cell or cluster-level scPrinter and seq2PRINT model training are not covered by this Skill's scripts.

Provide the BAMs per condition, consensus peaks, genome FASTA and build, blacklist, motif database, organism, and whether the data are bulk or single-cell. Ask for clarification when these choices change the tool.
