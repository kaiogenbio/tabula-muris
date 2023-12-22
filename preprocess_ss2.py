#!/usr/bin/env python

"""
preprocess_ss2.py

Preprocesses and generates reports on Smartseq2 (cell, gene_UMI_count) matrix.
"""

import argparse
import sys
import json

from collections import OrderedDict
from pathlib import Path
from io import BytesIO
from PIL import Image
from scipy.sparse import csr_matrix

from matplotlib.axes import Axes

import numpy as np
import scanpy as sc

CONSTANTS_SS2 = {
    "min_genes": 500, # taken from the Tabula Muris 2018 paper
    "min_reads": 5e4,
    "norm_count": 1e6
}


def sort_adata_by_obs_names(adata):
    return adata[np.argsort(adata.obs_names), :]


def preprocess(adata):
    # 1. QC metrics
    #    - # reads
    #    - # genes
    #    - % ERCC
    #    - [not available] % mitochondrial genes
    # 2. [not performed, optional] Correct for transcript length bias in SS2
    # 3. Normalize (reads per million)
    # 4. Log transform
    # 5. [not performed] Regress out technical covariates
    #    - e.g. total # of reads

    meta = {}
    axes = OrderedDict()

    print(f"Original:\n{adata}")

    # As outlined in the 2018 paper, we sort adata by adata.obs_names
    # for more reproducible results.
    adata = sort_adata_by_obs_names(adata)

    meta["orig_cell_names"] = sorted(adata.obs_names)
    meta["orig_gene_names"] = sorted(adata.var_names)

    # filter out low quality cells using cutoffs from the Tabula Muris 2018 paper
    sc.pp.filter_cells(adata, min_genes=CONSTANTS_SS2["min_genes"])
    print(f"Post min_genes:\n{adata}")

    sc.pp.filter_cells(adata, min_counts=CONSTANTS_SS2["min_reads"])
    print(f"Post min_reads:\n{adata}")

    # no mitochrondrial genes in this dataset
    # adata.var["ercc"] = adata.var_names.str.lower().str.startswith("ercc-")

    # instead we have ERCC spike-ins
    adata.var["ercc"] = adata.var_names.str.lower().str.startswith("ercc-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["ercc"], percent_top=None, log1p=False, inplace=True)

    # filter out ERCC genes
    adata = adata[:, ~adata.var["ercc"]]
    print(f"Post ERCC:\n{adata}")

    # QC plots
    axes["qc0"] = sc.pl.highest_expr_genes(
        adata,
        n_top=20,
        show=False)
    axes["qc1"] = sc.pl.violin(
        adata,
        "n_genes_by_counts",
        jitter=0.4,
        show=False)
    axes["total_counts"] = sc.pl.violin(
        adata,
        "total_counts",
        jitter=0.4,
        show=False)
    axes["ercc"] = sc.pl.violin(
        adata,
        "pct_counts_ercc",
        jitter=0.4,
        show=False)
    axes["qc2"] = sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_ercc",
        show=False)
    axes["qc3"] = sc.pl.scatter(adata,
        x="total_counts",
        y="n_genes_by_counts",
        show=False)

    # store a copy of counts data before normalization
    adata.layers["counts"] = csr_matrix(adata.X)

    # log-normalize counts
    sc.pp.normalize_total(adata, target_sum=CONSTANTS_SS2["norm_count"])
    sc.pp.log1p(adata)

    # store log-normalized counts in adata.raw
    adata.raw = adata

    meta["final_cell_names"] = sorted(adata.obs_names)
    meta["final_gene_names"] = sorted(adata.var_names)

    return adata, meta, axes


def get_prefix_from_ss2_csv(csv: Path):
    return csv.name.split("-")[0]


def get_all_ss2_csvs(csv_dir: Path):
    return csv_dir.glob("*.csv")


def convert_rgba_to_rgb(img: Image):
    if img.mode == "RGBA":
        rgb = Image.new("RGB", img.size, (255,255,255)) 
        rgb.paste(img, mask=img.split()[3])
        return rgb
    else:
        return img


def save_figs(axes: list[Axes], out_fname: Path):
    figs = []
    for ax in axes.values():
        buf = BytesIO()
        ax.get_figure().savefig(buf, transparent=False, dpi=600)
        buf.seek(0)
        img = Image.open(buf)
        img = convert_rgba_to_rgb(img)
        figs.append(img)

    figs[0].save(
        out_fname,
        "PDF",
        dpi=(600, 600),
        save_all=True,
        append_images=figs[1:]
    )

def main(*, outdir: Path, ss2_csv_dir: Path):
    stdout = sys.stdout
    sc.settings.verbosity = 3

    outdir.mkdir(parents=True, exist_ok=True)

    all_ss2_csv = get_all_ss2_csvs(ss2_csv_dir)

    for ss2_csv in all_ss2_csv:
        prefix = get_prefix_from_ss2_csv(ss2_csv)
        subdir = outdir / prefix
        subdir.mkdir(exist_ok=True)

        adata = sc.read_csv(ss2_csv, first_column_names=True).T

        with open(subdir / "preprocess.log", "w") as log:
            sys.stdout = log
            adata, meta, axes = preprocess(adata)
            sys.stdout = stdout
        
        with open(subdir / "metadata.json", "w") as fout:
            json.dump(meta, fout, indent=2)
        
        save_figs(axes, subdir / "figures.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--outdir",
        help="parent directory that will contain the subfolders with preprocessed SS2 data",
        required=True,
        type=Path
        )

    parser.add_argument(
        "--ss2dir",
        help="directory that contains all the Smartseq2 count CSV files",
        required=True,
        type=Path
        )

    args = parser.parse_args()

    main(outdir=args.outdir, ss2_csv_dir=args.ss2dir)
