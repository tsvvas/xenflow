#!/usr/bin/env python
import argparse
from pathlib import Path

from spatialdata_io import xenium


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--xenium_dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    sdata = xenium(args.xenium_dir.resolve())
    sdata.write(args.out.resolve())


if __name__ == "__main__":
    main()
