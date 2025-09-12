#!/usr/bin/env python
import argparse
from pathlib import Path

import geopandas as gpd
import scanpy as sc
import spatialdata as sd
from shapely.geometry import Polygon


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Split a SpatialData dataset by ROI polygons and write one .h5ad file per region."
    )
    p.add_argument(
        "--dataset-zarr",
        required=True,
        type=Path,
        help="SpatialData Zarr store that contains the expression table.",
    )
    p.add_argument(
        "--regions-zarr",
        required=True,
        type=Path,
        help="SpatialData Zarr store with `shapes['regions']` polygons.",
    )
    p.add_argument(
        "--sample-id",
        required=True,
        help="Primary sample ID (prepended to region_id in the output).",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    sdata = sd.read_zarr(args.dataset_zarr)
    sdata_regions = sd.read_zarr(args.regions_zarr)

    sdata.shapes["regions"] = sdata_regions.shapes["regions"]

    regions_gdf: gpd.GeoDataFrame = sd.transform(sdata.shapes["regions"], to_coordinate_system="global")

    for _, row in regions_gdf.iterrows():
        region_id = str(row["region_id"])
        polygon: Polygon = row.geometry

        roi_sdata = sdata.query.polygon(
            polygon=polygon,
            target_coordinate_system="global",
            filter_table=True,
            clip=False,
        )

        adata = roi_sdata.tables["xenium_table"]
        sc.pp.filter_cells(adata, min_counts=10)
        adata.obs["sample_uid"] = f"{args.sample_id}_{region_id}"

        out_path = Path(f"{args.sample_id}_{region_id}.h5ad").resolve()
        adata.write_h5ad(out_path)


if __name__ == "__main__":
    main()
