#!/usr/bin/env python
import argparse
import json
from functools import partial
from pathlib import Path

import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
import sopa
import sopa.aggregation.transcripts as tr
import spatialdata
from anndata import AnnData
from dask.diagnostics import ProgressBar
from scipy.sparse import coo_matrix, csr_matrix
from sopa import settings

sopa.settings.parallelization_backend = "dask"


# # Next two functions patch sopa to work with anndata >= 0.12.0
# def _add_coo(
#     X_partitions: list[coo_matrix],
#     geo_df: gpd.GeoDataFrame,
#     partition: pd.DataFrame,
#     gene_column: str,
#     gene_names: list[str],
# ) -> None:
#     if settings.gene_exclude_pattern is not None:
#         partition = partition[~partition[gene_column].str.match(settings.gene_exclude_pattern, case=False, na=False)]

#     points_gdf = gpd.GeoDataFrame(partition, geometry=gpd.points_from_xy(partition["x"], partition["y"]))
#     joined = geo_df.sjoin(points_gdf)
#     cells_indices, column_indices = joined.index, joined[gene_column].cat.codes

#     cells_indices = cells_indices[column_indices >= 0]
#     column_indices = column_indices[column_indices >= 0]

#     X_partition = coo_matrix(
#         (np.full(len(cells_indices), 1), (cells_indices, column_indices)),
#         shape=(len(geo_df), len(gene_names)),
#     )

#     X_partitions.append(X_partition)


# def _count_transcripts_aligned_csr(geo_df: gpd.GeoDataFrame, points: dd.DataFrame, value_key: str) -> AnnData:
#     """Count transcripts per cell. The cells and points have to be aligned (i.e., in the same coordinate system)

#     Args:
#         geo_df: Cells geometries
#         points: Transcripts dataframe
#         value_key: Key of `points` containing the genes names

#     Returns:
#         An `AnnData` object of shape `(n_cells, n_genes)` with the counts per cell
#     """
#     points[value_key] = points[value_key].astype("category").cat.as_known()
#     gene_names = points[value_key].cat.categories.astype(str)

#     X = csr_matrix((len(geo_df), len(gene_names)), dtype=int)
#     adata = AnnData(X=X, var=pd.DataFrame(index=gene_names))
#     adata.obs_names = geo_df.index.astype(str)

#     geo_df = geo_df.reset_index()

#     X_partitions = []

#     with ProgressBar():
#         points.map_partitions(
#             partial(_add_coo, X_partitions, geo_df, gene_column=value_key, gene_names=gene_names),
#             meta=(),
#         ).compute()

#     for X_partition in X_partitions:
#         adata.X += X_partition

#     if settings.gene_exclude_pattern is not None:
#         adata = adata[:, ~adata.var_names.str.match(settings.gene_exclude_pattern, case=False, na=False)].copy()

#     return adata


# tr._count_transcripts_aligned = _count_transcripts_aligned_csr


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="resegment_cells.py",
        description="Re-segment cells in a SpatialData Zarr using transcript-aware methods (Baysor/ProSeg) via SoPA.",
    )
    p.add_argument("--xenium_file", required=True, type=Path, help="Input SpatialData Zarr produced from Xenium.")
    p.add_argument("--out", required=True, type=Path, help="Output Zarr path to write (e.g., <sample>_reseg.zarr).")
    p.add_argument("--method", required=True, choices=["baysor", "proseg"], help="Segmentation method to use.")
    p.add_argument("--config", required=False, type=Path, help="Config file.")
    return p.parse_args()


def main():
    args = _parse_args()
    sdata = spatialdata.read_zarr(args.xenium_file.resolve())
    segment_config = json.loads(args.config.read_bytes())
    patches_config = {
        "patch_width": 2000,
        "patch_overlap": 50,
        "min_points_per_patch": 8000,
        "prior_shapes_key": "cell_boundaries",
    }

    if args.method == "baysor":
        sopa.make_transcript_patches(sdata, **patches_config)
        sopa.segmentation.baysor(sdata, config=segment_config, key_added=args.method, min_area=20)
    elif args.method == "proseg":
        sopa.make_transcript_patches(sdata, patch_width=None, prior_shapes_key="cell_boundaries")
        sopa.segmentation.proseg(sdata, key_added=args.method)

    sopa.aggregate(sdata, shapes_key=args.method, aggregate_channels=False, key_added=None)
    sdata.write(args.out.resolve())


if __name__ == "__main__":
    main()
