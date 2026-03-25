#!/usr/bin/env python
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg

RAW_COUNTS = "counts"
LOGSCALED_COUNTS = "scaled_counts"
CLUSTER_COL = "leiden"
UNASSIGNED_KEY = "Unassigned"


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


def preprocess_standard(adata):
    sc.pp.filter_cells(adata, min_counts=10)
    adata.layers[RAW_COUNTS] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers[LOGSCALED_COUNTS] = adata.X.copy()
    sc.pp.highly_variable_genes(adata)
    return adata


def cluster_adata(adata):
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor="igraph")
    return adata


def get_marker_genes(adata_in, cluster_col, top_n):
    adata = adata_in.copy()
    adata.X = adata.layers[LOGSCALED_COUNTS]
    sc.tl.rank_genes_groups(adata, groupby=cluster_col, use_raw=False, method="wilcoxon")
    groups = adata.obs[cluster_col].astype("category").cat.categories
    markers = []
    for group in groups:
        df = sc.get.rank_genes_groups_df(adata, group=group)
        top_genes = df.sort_values("scores", ascending=False)["names"].head(top_n).tolist()
        markers.extend(top_genes)
    return markers


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
    adata_st = adata_xm.copy()
    adata_st = preprocess_standard(adata_st)
    adata_st = cluster_adata(adata_st)

    adata_sc = sc.read(args.reference_file.resolve())
    adata_sc = preprocess_standard(adata_sc)

    genes_sc = get_marker_genes(adata_sc, args.cell_type_key, 10)
    genes_st = get_marker_genes(adata_st, CLUSTER_COL, 10)

    genes_shared = np.intersect1d(genes_sc, genes_st)

    if genes_shared.shape[0] == 0:
        raise ValueError("0 genes selected for mapping. Adjust the config or find a better reference.")

    tg.pp_adatas(adata_sc, adata_st, genes=genes_shared, gene_to_lowercase=False)

    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_st,
        mode="clusters",
        cluster_label=args.cell_type_key,
        density_prior="uniform",
        num_epochs=500,
        device="cuda",
    )

    tg.plot_training_scores(ad_map, bins=20, alpha=0.5)
    plt.savefig(args.training_plot, dpi=300, bbox_inches="tight")
    plt.close()

    tg.project_cell_annotations(ad_map, adata_st, annotation=args.cell_type_key)
    adata_st = get_tangram_predictions(adata_st)

    tangram_df = adata_xm.obs.join(adata_st.obs[["tangram_prediction", "tangram_score"]], how="left")
    adata_xm.obs = tangram_df
    adata_xm.obs["tangram_prediction"] = adata_xm.obs["tangram_prediction"].fillna(UNASSIGNED_KEY)

    adata_st.write_h5ad(args.out.resolve())


if __name__ == "__main__":
    main()
