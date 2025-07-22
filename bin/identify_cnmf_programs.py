#!/usr/bin/env python
"""
identify_programs.py
Run consensus-NMF on one .h5ad file and keep **all** cNMF output.

Output layout
-------------
<out_dir>/                       (passed via --out-dir)
    raw/ …                       # all intermediate runs
    k_selection.png              # elbow-plot
    gene_spectra_tpm.k*.csv      # one per k
    cell_usage.k*.csv            # one per k
"""

import argparse
import re
from pathlib import Path

import numpy as np
from cnmf import cNMF


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, type=Path, help=".h5ad produced by SPLIT_SAMPLES")
    ap.add_argument("--out-dir", required=True, type=Path, help="Directory to create (will be <uid>_cnmf)")
    ap.add_argument(
        "--components", default="10,15,20", type=parse_comma_ints, help="Comma-separated list of k; e.g. 5,10,15"
    )
    ap.add_argument("--n-iter", type=int, default=100)
    ap.add_argument("--num-highvar-genes", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=4)
    return ap.parse_args()


def parse_comma_ints(text: str) -> tuple[int, ...]:
    """
    Convert '5, 15, 5'  →  (5, 15, 5)
    Converts any amount of whitespace and validates that every
    chunk is an integer.
    """
    parts = re.split(r"\s*,\s*", text.strip())
    try:
        return tuple(map(int, parts))
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"components must be ints: {text!r}") from e


def main() -> None:
    args = _parse_args()
    k_list = np.arange(*args.components)

    cnmf = cNMF(name=args.out_dir.as_posix())
    cnmf.prepare(
        counts_fn=args.h5ad.as_posix(),
        components=k_list,
        n_iter=args.n_iter,
        num_highvar_genes=args.num_highvar_genes,
        seed=args.seed,
    )
    cnmf.factorize_multi_process(total_workers=args.workers)
    cnmf.combine()
    cnmf.k_selection_plot()


if __name__ == "__main__":
    main()
