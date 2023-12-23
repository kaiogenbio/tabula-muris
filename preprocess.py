#!/usr/bin/env python

from collections import OrderedDict
import scanpy as sc

from scipy.sparse import csr_matrix
from jsonschema import validate

import utils


def preprocess_adata(adata, constants):
    constants_schema = {
        "type": "object",
        "properties": {
            "min_genes": {"type": "number"},
            "min_counts": {"type": "number"}, # for UMI based datasets (i.e. 10X, min_counts is actually min_UMI [UMIs collapse down PCR amplification bias whereas])
            "norm_total": {"type": "number"}
        },
        "required": ["min_genes", "min_counts", "norm_total"],
        "additionalProperties": False
    }

    validate(instance=constants, schema=constants_schema)

    meta = {}
    axes = OrderedDict()

    print(f"Original: {adata}")

    meta["orig_cell_names"] = sorted(adata.obs_names)
    meta["orig_gene_names"] = sorted(adata.var_names)

    adata = utils.sort_adata_by_obs_names(adata)

    # 1. filter out low quality cells based on low UMIs and low #genes
    sc.pp.filter_cells(adata, min_genes=constants["min_genes"])
    print(f"Post min_genes: {adata}")

    sc.pp.filter_cells(adata, min_counts=constants["min_counts"])
    print(f"Post min_counts/umi: {adata}")

    # 2. use QC filters to filter out cells with high # of ERCC reads
    # NB: 10X sequencing doesn't appear to have any ERCC spike-ins
    adata.var["ercc"] = adata.var_names.str.lower().str.startswith("ercc-")
    print("Num ERCC genes:", sum(adata.var["ercc"]))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["ercc"],
        percent_top=None,
        log1p=False,
        inplace=True)

    # filter out ERCC genes
    adata = adata[:, ~adata.var["ercc"]]

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

    adata.layers["counts"] = csr_matrix(adata.X.copy())

    # 3. normalize total, take log
    sc.pp.normalize_total(adata, target_sum=constants["norm_total"])
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = csr_matrix(adata.X.copy())

    adata.X = csr_matrix(adata.X)

    adata.raw = adata

    meta["final_cell_names"] = sorted(adata.obs_names)
    meta["final_num_cells"] = len(adata.obs_names)
    meta["final_gene_names"] = sorted(adata.var_names)
    meta["final_num_genes"] = len(adata.var_names)

    return adata, meta, axes
