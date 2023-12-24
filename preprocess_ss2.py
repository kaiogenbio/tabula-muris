#!/usr/bin/env python

"""
preprocess_ss2.py

Preprocesses and generates reports on Smartseq2 (cell, gene_UMI_count) matrix.
"""

import argparse
import json

from pathlib import Path

import scanpy as sc
import pandas as pd

from preprocess import preprocess_adata
import utils


# constants are taken from Tabula Muris 2018 paper
CONSTANTS_SS2 = {
    "min_genes": 500,
    "min_counts": 5e4,
    "norm_total": 1e6
}


def get_prefix_from_ss2_csv(csv: Path):
    return csv.name.rsplit("-", 1)[0]


def get_all_ss2_csvs(csv_dir: Path):
    return csv_dir.glob("*.csv")


def main(*, outdir: Path, ss2_csv_dir: Path, annotations_file: Path):
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
            "Neurog3>0_raw",
            "Neurog3>0_scaled",
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
            "subsetE",
            "subsetE_cluster.ids",
            "subtissue"
            ],
        inplace=True
        )
    sc.settings.verbosity = 3

    all_ss2_csv = get_all_ss2_csvs(ss2_csv_dir)

    for ss2_csv in all_ss2_csv:
        prefix = get_prefix_from_ss2_csv(ss2_csv)
        subdir = outdir / prefix
        subdir.mkdir(parents=True, exist_ok=True)

        print(f"Processing {prefix}...")

        adata = sc.read_csv(ss2_csv, first_column_names=True).T
        adata.obs_names_make_unique()

        with utils.Tee([open(subdir / "preprocess.log", "w")]):
            adata, meta, axes = preprocess_adata(
                adata, CONSTANTS_SS2)

        # merge annotations to adata.obs
        adata.obs = adata.obs.join(annotations)

        # remove cells without a cell_ontology_class assigned
        adata = adata[~adata.obs["cell_ontology_class"].isna(), :]
        
        with open(subdir / "metadata.json", "w") as fout:
            json.dump(meta, fout, indent=2)
        
        utils.save_axes(axes, subdir / "figures.pdf")
        adata.write_h5ad(subdir / "preprocessed.h5ad")


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

    parser.add_argument(
        "--annotations",
        help="annotations file containing expert cell ontology labels",
        required=True,
        type=Path
    )

    args = parser.parse_args()

    main(
        outdir=args.outdir,
        ss2_csv_dir=args.ss2dir,
        annotations_file=args.annotations
        )
