#!/usr/bin/env python

import argparse
import json

from pathlib import Path

import pandas as pd
import scanpy as sc

from preprocess import preprocess_adata
import utils


CONSTANTS_10X = {
    "min_genes": 500,
    "min_counts": 1000,
    "norm_total": 1e4
}


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
            adata, meta, axes = preprocess_adata(
                adata, CONSTANTS_10X
                )

        # save original indices as "barcode" column
        adata.obs["barcode"] = adata.obs.index
        # match index with annotations_file for join
        utils.modify_index(adata.obs, lambda x: f"{channel}_{x[:-2]}")
        # merge annotations to adata.obs
        adata.obs = adata.obs.join(annotations)

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
