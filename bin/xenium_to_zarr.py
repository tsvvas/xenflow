#!/usr/bin/env python
import argparse
from pathlib import Path

from spatialdata_io import xenium

DEFAULT_XENFLOW_CONFIG = {
    "current": {
        "morphology_image": "morphology_focus",
        "nucleus_shapes": "nucleus_boundaries",
        "cell_shapes": "cell_boundaries",
        "cell_table": "xenium_table",
    },
    "history": [],
    "schema_version": "",
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--xenium_dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    sdata = xenium(args.xenium_dir.resolve())
    sdata.tables["xenium_table"] = sdata.tables["table"].copy()
    sdata.attrs["xenflow"] = DEFAULT_XENFLOW_CONFIG
    sdata.write(args.out.resolve())


if __name__ == "__main__":
    main()
