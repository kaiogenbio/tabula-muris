#!/usr/bin/env python

import argparse
import sys
import json

from pathlib import Path
from collections import OrderedDict

import pandas as pd
import numpy as np
import scanpy as sc

from scipy.sparse import csr_matrix

import utils

CONSTANTS_10X = {
    "min_genes": 501,
    "min_umi": 1001,
    "norm_total": 1e4
}

def preprocess(adata):
    meta = {}
    axes = OrderedDict()

    adata.layers["counts"] = csr_matrix(adata.X.copy())

    print("Original adata:", adata)

    # 1. filter out low quality cells based on low UMIs and low #genes
    sc.pp.filter_cells(adata, min_genes=CONSTANTS_10X["min_genes"])
    sc.pp.filter_cells(adata, min_counts=CONSTANTS_10X["min_umi"])

    print("After initial filtering steps:", adata)

    # 2. use QC filters to filter out cells with high # of ERCC reads
    # NB: 10X sequencing doesn't appear to have any ERCC spike-ins
    adata.var["ercc"] = adata.var_names.str.lower().str.startswith("ercc-")
    print("Num ERCC genes:", sum(adata.var["ercc"]))
    adata = adata[:, ~adata.var["ercc"]]

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["ercc"],
        percent_top=None,
        log1p=False,
        inplace=True)
    axes["highest expressed genes"] = sc.pl.highest_expr_genes(
        adata,
        n_top=20,
        show=False)
    axes["num unique genes expressed"] = sc.pl.violin(
        adata,
        "n_genes_by_counts",
        jitter=0.4,
        show=False)
    axes["total_counts"] = sc.pl.violin(
        adata,
        "total_counts",
        jitter=0.4,
        show=False)
    axes["pct_counts_ercc"] = sc.pl.violin(
        adata,
        "pct_counts_ercc",
        jitter=0.4,
        show=False)
    axes["total_counts vs. pct_counts_ercc"] = sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_ercc",
        show=False)
    axes["total_counts vs. n_genes_by_counts"] = sc.pl.scatter(adata,
        x="total_counts",
        y="n_genes_by_counts",
        show=False)

    # 3. normalize total, take log
    sc.pp.normalize_total(adata, target_sum=CONSTANTS_10X["norm_total"])
    sc.pp.log1p(adata)

    adata.X = csr_matrix(adata.X)
    adata.raw = adata

    meta["cell_names"] = sorted(adata.obs_names)
    meta["num_cells"] = len(adata.obs_names)
    meta["gene_names"] = sorted(adata.var_names)
    meta["num_genes"] = len(adata.var_names)

    return adata, meta, axes


def main(*, outdir, dropletdir, annotations_file: Path):
    annotations = pd.read_csv(
        annotations_file,
        index_col="cell",
        dtype={
            "cell": "object",
            "cell_ontology_class": "object",
            "cell_ontology_id": "object",
            "channel": "object",
            "cluster.ids": "object",
            "mouse.id": "object",
            "mouse.sex": "object",
            "subsetA": "object",
            "subsetA_cluster.ids": "object",
            "subsetB": "object",
            "subsetB_cluster.ids": "object",
            "subsetC": "object",
            "subsetC_cluster.ids": "object",
            "subsetD": "object",
            "subsetD_cluster.ids": "object",
            "subtissue": "object",
            "tissue": "object"
        })
    annotations.drop(
        columns=[
            "tissue_tSNE_1",
            "tissue_tSNE_2",
            "free_annotation",
            "subsetA",
            "subsetA_cluster.ids",
            "subsetB",
            "subsetB_cluster.ids",
            "subsetC",
            "subsetC_cluster.ids",
            "subsetD",
            "subsetD_cluster.ids",
            "subtissue"
            ],
        inplace=True
        )

    assert annotations.index.is_unique, f"{annotations_file} index using index_col='cell' is not unique."

    all_subdirs = utils.get_subdirs(dropletdir)

    for subdir in all_subdirs:
        prefix = subdir.name

        out_subdir = outdir / prefix
        out_subdir.mkdir(parents=True, exist_ok=True)

        print(f"Processing {prefix}...")

        tissue, channel = prefix.split("-")
        adata = sc.read_10x_mtx(subdir)
        adata.obs_names_make_unique()

        # perform QC
        with utils.Tee([open(out_subdir / "preprocess.log", "w")]):
            adata, meta, axes = preprocess(adata)

        # save original indices as "barcode" column
        adata.obs["barcode"] = adata.obs.index
        # match index with annotations_file for join
        utils.modify_index(adata.obs, lambda x: f"{channel}_{x[:-2]}")
        # merge annotations to adata.obs
        adata.obs = adata.obs.join(annotations, validate="1:1")

        # remove cells without a cell_ontology_class assigned
        adata = adata[~adata.obs["cell_ontology_class"].isna(), :]

        with open(out_subdir / "metadata.json", "w") as fout:
            json.dump(meta, fout, indent=2)

        utils.save_axes(axes, out_subdir / "figures.pdf")
        adata.write_h5ad(out_subdir / "preprocessed.h5ad")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--dropletdir", required=True, type=Path)
    parser.add_argument("--annotations", required=True, type=Path)

    args = parser.parse_args()

    main(
        outdir=args.outdir,
        dropletdir=args.dropletdir,
        annotations_file=args.annotations
        )
