#!/usr/bin/env python

import streamlit as st
import scanpy as sc
import pandas as pd

def prep_adata(adata, rep=None):
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30 if rep is None else None, use_rep=rep)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)


def show_plots(adata):
    sc.set_figure_params(dpi=150, figsize=(8, 8))
    sc.pl.umap(adata, color=["leiden"], frameon=False)
    sc.pl.umap(adata, color=["tissue"], frameon=False)
    sc.pl.umap(adata, color=["channel"], frameon=False)
    sc.pl.umap(adata, color=["cell_ontology_class"], frameon=False, legend_loc="right margin", legend_fontsize="xx-small", legend_fontweight="light")


def plot_adata(adata, obsm_key, color_key):
    df = pd.DataFrame({
        "umap1": adata.obsm[obsm_key][:, 0],
        "umap2": adata.obsm[obsm_key][:, 1],
        "color": adata.obs[color_key]
        }, index=adata.obs.index)
    
    st.scatter_chart(
        df,
        x="umap1",
        y="umap2",
        size=10,
        color="color")


def main():
    adata_concat = sc.read_h5ad("../integrated/10x/concatenated.h5ad")
    adata_scvi = sc.read_h5ad("../integrated/10x/scvi-integrated.h5ad")

    prep_adata(adata_concat)
    prep_adata(adata_scvi)

    plot_adata(adata_scvi, "X_umap", "cell_ontology_class")


if __name__ == "__main__":
    main()