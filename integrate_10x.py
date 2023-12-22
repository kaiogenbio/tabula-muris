#!/usr/bin/env python

from typing import List
from pathlib import Path

import argparse

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


def main(*, glob_str, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)

    adata_paths = sorted(Path(".").glob(glob_str))
    all_adata = concat_adatas(adata_paths)
    sparsify(all_adata)

    all_adata.write_h5ad(out_dir / "integrated.h5ad")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--glob", required=True)
    parser.add_argument("--outdir", required=True, type=Path)

    args = parser.parse_args()

    main(glob_str=args.glob, out_dir=args.outdir)