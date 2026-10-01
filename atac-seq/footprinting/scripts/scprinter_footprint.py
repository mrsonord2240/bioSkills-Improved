#!/usr/bin/env python3
"""Bulk multi-scale footprint scores with scPrinter 1.2.0 (classic get_footprint_score route).

Steps: build a custom Genome, predict the genome-wide Tn5 bias (GPU strongly advised), import one fragment file as a single
pseudobulk sample, then score every region in --regions at each footprint scale (mode).

Usage:
  scprinter_footprint.py --fragments frags.tsv.gz --fasta hg38.fa --gtf genes.gtf --blacklist bl.bed \
      --regions sites.bed --outdir out [--bias bias.h5] [--modes 2-100] [--width 200] [--device cuda:0]

--fragments  bgzip + tabix BED with columns chrom, start, end, barcode; raw (unshifted) fragment ends,
             e.g. from properly paired reads: chrom, POS-1, POS-1+TLEN, name. scPrinter applies the +4/-5 Tn5 shift.
--regions    BED; each region is resized to --width around its centre (motif centres for footprint work).
Outputs: <outdir>/footprints.npz (scores regions x modes x width, keys, modes) and <outdir>/center_by_mode.tsv.
Set SCPRINTER_DATA to a writable directory: the first import downloads the pretrained models (needs internet).
"""
import argparse
import os
import sys


def parse_modes(text):
    import numpy as np
    if "-" in text:
        lo, hi = (int(x) for x in text.split("-"))
        return np.arange(lo, hi + 1)
    return np.array([int(x) for x in text.split(",")])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--blacklist", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bias", help="existing genome Tn5 bias .h5; predicted into <outdir>/bias.h5 when absent")
    ap.add_argument("--modes", default="2-100", help="scales: 'lo-hi' or comma list (default 2-100)")
    ap.add_argument("--width", type=int, default=200)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--name", default="sample")
    ap.add_argument("--min-fragments", type=int, default=1000)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    os.environ.setdefault("SCPRINTER_DATA", os.path.join(os.path.abspath(a.outdir), "scprinter_data"))
    os.makedirs(os.environ["SCPRINTER_DATA"], exist_ok=True)
    import numpy as np
    import pandas as pd
    import scprinter as scp

    modes = parse_modes(a.modes)
    bias = a.bias or os.path.join(a.outdir, "bias.h5")
    if not os.path.exists(bias):
        scp.genome.predict_genome_tn5_bias(fa_file=a.fasta, save_name=bias, tn5_model=scp.datasets.pretrained_Tn5_bias_model,
                                           context_radius=50, device=a.device, batch_size=5000)
    genome = scp.genome.Genome(name="custom", fa_file=a.fasta, gff_file=a.gtf, bias_file=bias, blacklist_file=a.blacklist)
    h5 = os.path.join(a.outdir, "printer.h5ad")
    if os.path.exists(h5):
        os.remove(h5)
    printer = scp.pp.import_fragments(path_to_frags=[a.fragments], barcodes=[None], savename=h5, genome=genome,
                                      sample_names=[a.name], min_num_fragments=a.min_fragments, min_tsse=0,
                                      sorted_by_barcode=False, low_memory=False)
    printer.load_disp_model()
    regions = pd.read_csv(a.regions, sep="\t", header=None, comment="#").iloc[:, :3]
    regions.columns = ["chrom", "start", "end"]
    scp.tl.get_footprint_score(printer, [list(printer.obs_names)], [a.name], regions, region_width=a.width, modes=modes,
                               footprintRadius=None, flankRadius=None, n_jobs=a.jobs, save_key="fp", backed=True,
                               overwrite=True)
    ad = printer.footprintsadata["fp"]
    keys = list(ad.obsm.keys())                           # "chrom:start-end" of the resized regions, same order as scores
    scores = np.stack([np.asarray(ad.obsm[k])[0] for k in keys])
    if not np.isfinite(scores).all():
        sys.exit("ERROR: non-finite footprint scores")
    np.savez_compressed(os.path.join(a.outdir, "footprints.npz"), scores=scores, keys=np.array(keys), modes=modes)
    mid = a.width // 2
    center = scores[:, :, mid - 10:mid + 10].mean(axis=2).mean(axis=0)
    flank = np.concatenate([scores[:, :, :mid - 40], scores[:, :, mid + 40:]], axis=2).mean(axis=2).mean(axis=0)
    pd.DataFrame({"mode": modes, "center_pm10bp_mean": center, "flank_mean": flank}).to_csv(
        os.path.join(a.outdir, "center_by_mode.tsv"), sep="\t", index=False)
    print(f"scprinter {scp.__version__}; scores {scores.shape} (regions x modes x bp); wrote {a.outdir}/footprints.npz")


if __name__ == "__main__":
    main()
