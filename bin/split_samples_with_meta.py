#!/usr/bin/env python
import argparse
import json
import warnings
from pathlib import Path

import spatialdata as sd
from shapely.geometry import box


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split a SpatialData dataset into one .h5ad per region.")
    parser.add_argument("--xenium-file", required=True, type=Path, help="SpatialData Zarr store.")
    parser.add_argument("--sample-id", required=True, help="Sample/slide ID to look up regions in the metadata JSON.")
    parser.add_argument(
        "--metadata-file",
        required=True,
        type=Path,
        help="JSON metadata file containing sample ID -> region bbox mappings.",
    )
    return parser.parse_args()


def get_region_data(meta, sample_id):
    matches = [k for k, v in meta.items() if v["slide_id"] == sample_id]
    if not matches:
        raise ValueError(f"No metadata entry found for sample_id={sample_id!r}")
    if len(matches) > 1:
        raise ValueError(f"Multiple metadata entries found for sample_id={sample_id!r}: {matches}")
    record = matches[0]
    return meta[record]["sample_coords"]


def main():
    args = _parse_args()
    sdata = sd.read_zarr(args.xenium_file)
    current_table = sdata.attrs["xenflow"]["current"]["tx_table"]
    current_nucleus = sdata.attrs["xenflow"]["current"]["nucleus_shapes"]
    sdata_filtered = sdata.subset(element_names=[current_nucleus, current_table])

    meta = json.loads(args.metadata_file.read_text(encoding="utf-8"))
    regions = get_region_data(meta, args.sample_id)

    for region_name, region_data in regions.items():
        new_sample_id = f"{args.sample_id}_{region_name}"
        bbox = box(*region_data["coordinates"])

        roi_sdata = sdata_filtered.query.polygon(
            polygon=bbox,
            target_coordinate_system="global",
            filter_table=True,
            clip=False,
        )

        adata = roi_sdata.tables[current_table].copy()
        adata.obs["sample_uid"] = new_sample_id
        adata.obs["patient_id"] = region_data.get("patient_id")
        adata.obs["tissue_id"] = region_data.get("tissue_id")

        out_path = Path(f"{new_sample_id}.h5ad").resolve()
        warnings.warn(f"{out_path!r}: {adata.n_obs} cells, {adata.n_vars} genes.")
        adata.write_h5ad(out_path)


if __name__ == "__main__":
    main()
