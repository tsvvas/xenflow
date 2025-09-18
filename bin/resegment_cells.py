#!/usr/bin/env python
import argparse
import json
import os
from pathlib import Path

import sopa
import spatialdata

sopa.settings.parallelization_backend = "dask"
sopa.settings.dask_client_kwargs = {
    "n_workers": int(os.getenv("REQUESTED_CPUS")),
    "processes": True,
    "threads_per_worker": 1,
    "dashboard_address": ":8787",
}
sopa.log.setLevel(sopa.logging.ERROR)


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
        "patch_width": 1000,
        "patch_overlap": 50,
        "min_points_per_patch": 8000,
        "prior_shapes_key": "cell_boundaries",
    }

    sdata.tables["xenium_table"] = sdata.tables["table"].copy()

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
