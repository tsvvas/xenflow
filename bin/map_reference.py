#!/usr/bin/env python
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg

LAYER_KEY = "counts"


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--xenium-file", required=True, type=Path)
    p.add_argument("--reference-file", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument(
        "--cell-type-key",
        required=True,
        type=str,
        help="Column in reference AnnData .obs containing cell type labels.",
    )
    p.add_argument(
        "--training-plot",
        default=None,
        required=True,
        type=Path,
        help="Optional path to save Tangram training scores plot.",
    )
    args = p.parse_args()
    return args


def cluster_adata(adata, hvg=True, n_top_genes=2000):
    sc.pp.filter_cells(adata, min_counts=10)
    adata.layers[LAYER_KEY] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    if hvg:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor="igraph")
    return adata


def bin_entropy(x, axis=1):
    """Calculates Shannon entropy with the base of 2 on count data."""
    x = np.atleast_2d(x).astype(int)

    if axis == 0:
        x = x.T

    nrows, ncols = x.shape
    nbins = x.max() + 1

    counts = np.vstack([np.bincount(row, minlength=nbins) for row in x])

    p = counts / float(ncols)
    p = np.where(p > 0, p, 1)

    return -np.sum(p * np.log2(p), axis=1)


def calculate_cluster_entropy(adata, cluster_col):
    X = adata.X.toarray()
    cluster_labels = adata.obs[cluster_col].values
    unique_clusters = np.unique(cluster_labels)
    n_clusters = len(unique_clusters)

    avg_gene_values = np.zeros((n_clusters, X.shape[1]))

    for i, cluster in enumerate(unique_clusters):
        cluster_indices = np.where(cluster_labels == cluster)[0]

        avg_gene_values[i, :] = X[cluster_indices].mean(axis=0)

    avg_gene_values_transposed = avg_gene_values.T

    entropies = bin_entropy(avg_gene_values_transposed, axis=1)

    return entropies


def get_informative_genes(adata, cluster_col, threshold=0.01):
    gene_names = adata.var_names
    entropies = calculate_cluster_entropy(adata, cluster_col)
    gene_entropies = zip(gene_names, entropies)
    filtered_genes = filter(lambda x: x[1] > threshold, gene_entropies)
    return [gene[0] for gene in filtered_genes]


def get_tangram_predictions(adata, perc=0, merge=True):
    tg_ct_pred = adata.obsm["tangram_ct_pred"]
    tg_ct_pred = tg_ct_pred.clip(tg_ct_pred.quantile(perc), tg_ct_pred.quantile(1 - perc), axis=1)
    tg_ct_pred = (tg_ct_pred - tg_ct_pred.min()) / (tg_ct_pred.max() - tg_ct_pred.min())

    df = pd.DataFrame()
    df["tangram_prediction"] = tg_ct_pred.idxmax(axis=1)
    df["tangram_score"] = tg_ct_pred.apply(lambda x: max(x), axis=1)

    if merge:
        adata.obs = adata.obs.join(df)
        return adata

    return df


def main():
    args = _parse_args()
    adata_xm = sc.read(args.xenium_file.resolve())
    adata_sc = sc.read(args.reference_file.resolve())

    adata = adata_xm.copy()
    adata = cluster_adata(adata)

    genes_sc = get_informative_genes(adata_sc, args.cell_type_key)
    genes_st = get_informative_genes(adata, "leiden")

    genes_shared = np.intersect1d(genes_sc, genes_st)
    tg.pp_adatas(adata_sc, adata, genes=genes_shared, gene_to_lowercase=False)

    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata,
        mode="clusters",
        cluster_label=args.cell_type_key,
        density_prior="uniform",
        num_epochs=500,
        device="cuda",
    )

    tg.plot_training_scores(ad_map, bins=20, alpha=0.5)
    plt.savefig(args.training_plot, dpi=300, bbox_inches="tight")
    plt.close()

    tg.project_cell_annotations(ad_map, adata, annotation=args.cell_type_key)
    adata = get_tangram_predictions(adata)

    tangram_df = adata_xm.obs.join(adata.obs[["tangram_prediction", "tangram_score"]], how="left")
    adata_xm.obs = tangram_df
    adata_xm.obs["tangram_prediction"] = adata_xm.obs["tangram_prediction"].fillna("unassigned")

    adata_xm.write_h5ad(args.out.resolve())


if __name__ == "__main__":
    main()
