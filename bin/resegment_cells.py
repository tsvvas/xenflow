#!/usr/bin/env python
import argparse
import json
from pathlib import Path

import geopandas as gpd
import sopa
import spatialdata


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="resegment_cells.py",
        description="Re-segment cells in a SpatialData Zarr using transcript-aware methods (Baysor/ProSeg) via SoPA.",
    )
    p.add_argument("--xenium_file", required=True, type=Path, help="Input SpatialData Zarr produced from Xenium.")
    p.add_argument(
        "--out", required=True, type=Path, help="Output Zarr path to write (e.g., <sample>_reseg_cells.zarr)."
    )
    p.add_argument("--method", required=True, choices=["baysor", "proseg"], help="Segmentation method to use.")
    p.add_argument(
        "--segment-config",
        dest="segment_config",
        required=True,
        type=Path,
        help="JSON file with method-specific segmentation parameters.",
    )
    p.add_argument(
        "--patch-config",
        dest="patch_config",
        required=True,
        type=Path,
        help="JSON file with transcript patching parameters.",
    )
    p.add_argument("--n-workers", default=1, type=int, required=False, help="Number of Dask workers to use.")
    return p.parse_args()


def setup_sopa(n_workers: int):
    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs = {
        "n_workers": n_workers,
        "processes": True,
        "threads_per_worker": 1,
        "dashboard_address": ":8787",
    }
    sopa.log.setLevel(sopa.logging.ERROR)


def filter_cells(
    sdata,
    nucleus_key="cellpose_boundaries",
    cell_key="baysor",
    coordinate_system="global",
):
    """Filters cell shapes that do not contain a nucleus centroid. Saves nuclei centroids in sdata.points and
    filtered shapes in sdata.shapes.
    """
    out_cells_key = cell_key + "_filtered"
    out_points_key = nucleus_key + "_centroids"

    nucs_gdf = spatialdata.transform(sdata.shapes[nucleus_key], to_coordinate_system=coordinate_system).copy()
    cells_gdf = spatialdata.transform(sdata.shapes[cell_key], to_coordinate_system=coordinate_system).copy()

    cent_df = spatialdata.get_centroids(nucs_gdf, coordinate_system=coordinate_system)
    cent_df = cent_df.compute() if hasattr(cent_df, "compute") else cent_df
    cent_df = cent_df.copy()
    cent_df["nucs_id"] = cent_df.index.astype(str)

    cent_gdf = gpd.GeoDataFrame(
        cent_df,
        geometry=gpd.points_from_xy(cent_df["x"], cent_df["y"]),
        index=cent_df.index,
    )

    cells_gdf = cells_gdf.copy()
    cells_gdf["cell_id"] = cells_gdf.index.astype(str)

    cent_in_cell = gpd.sjoin(
        cent_gdf[["nucs_id", "geometry"]],
        cells_gdf[["cell_id", "geometry"]],
        predicate="within",
        how="inner",
    )
    cell_ids_keep = cent_in_cell["cell_id"].unique()

    cells_filt_gdf = cells_gdf.loc[cells_gdf["cell_id"].isin(cell_ids_keep)].copy()
    cells_filt_gdf = cells_filt_gdf.set_index("cell_id", drop=True)

    cells_filt_elem = spatialdata.models.ShapesModel.parse(cells_filt_gdf)

    points_elem = spatialdata.models.PointsModel.parse(cent_df, coordinates={"x": "x", "y": "y"})

    spatialdata.transformations.set_transformation(
        cells_filt_elem, spatialdata.transformations.Identity(), coordinate_system
    )
    spatialdata.transformations.set_transformation(
        points_elem, spatialdata.transformations.Identity(), coordinate_system
    )

    sdata.shapes[out_cells_key] = cells_filt_elem
    sdata.points[out_points_key] = points_elem

    return sdata


def main():
    args = _parse_args()
    setup_sopa(args.n_workers)

    sdata = spatialdata.read_zarr(args.xenium_file.resolve())
    segment_config = json.loads(args.segment_config.read_bytes())
    patches_config = json.loads(args.patch_config.read_bytes())

    prior_shapes_key = sdata.attrs["xenflow"]["current"]["nucleus_shapes"]

    if args.method == "baysor":
        sopa.make_transcript_patches(sdata, prior_shapes_key=prior_shapes_key, **patches_config)
        sopa.segmentation.baysor(sdata, config=segment_config, key_added=args.method, min_area=20)
        sdata = filter_cells(sdata, nucleus_key=prior_shapes_key, cell_key=args.method)
    elif args.method == "proseg":
        sopa.make_transcript_patches(sdata, patch_width=None, prior_shapes_key=prior_shapes_key)
        sopa.segmentation.proseg(sdata, key_added=args.method)

    current_cell_shapes_key = args.method + "_filtered"
    sdata.attrs["xenflow"]["current"]["cell_shapes"] = current_cell_shapes_key

    image_key = sdata.attrs["xenflow"]["current"]["morphology_image"]
    sopa.aggregate(
        sdata, image_key=image_key, shapes_key=current_cell_shapes_key, aggregate_channels=False, key_added=None
    )
    sdata.attrs["xenflow"]["current"]["tx_table"] = current_cell_shapes_key + "_table"
    sdata.write(args.out.resolve())


if __name__ == "__main__":
    main()
