#!/usr/bin/env python

from typing import List
from pathlib import Path

import argparse

import scvi
import scanpy as sc

from scipy.sparse import csr_matrix


def concat_adatas(path_list: List[Path]):
    all_adata = None

    for adata_path in path_list:
        adata = sc.read_h5ad(adata_path)

        if all_adata is None:
            all_adata = adata
        else:
            all_adata = sc.concat([all_adata, adata])
    
    return all_adata


def sparsify(adata):
    for key, layer in adata.layers.items():
        adata.layers[key] = csr_matrix(layer)

    adata.X = csr_matrix(adata.X)


def main(*, glob_str, out_dir, batch_key):
    out_dir.mkdir(parents=True, exist_ok=True)

    adata_paths = sorted(Path(".").glob(glob_str))
    adata = concat_adatas(adata_paths)
    sparsify(adata)

    adata.write_h5ad(out_dir / "concatenated.h5ad")

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_key, # "channel", "plate.barcode"
        continuous_covariate_keys=None,
        categorical_covariate_keys=None)
    model = scvi.model.SCVI(adata)
    model.train()

    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.obsm["X_scVI_MDE"] = scvi.model.utils.mde(model.get_latent_representation())

    model.save(out_dir / "scvi-model.pt", overwrite=True)

    adata.write_h5ad(out_dir / "scvi-integrated.h5ad")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--glob", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--batchkey", required=True)

    args = parser.parse_args()

    main(glob_str=args.glob, out_dir=args.outdir, batch_key=args.batchkey)