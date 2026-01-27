#!/usr/bin/env python
import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import sopa
import spatialdata
from scipy.optimize import nnls
from skimage.exposure import rescale_intensity
from skimage.morphology import disk, white_tophat
from skimage.restoration import rolling_ball

DAPI_KEY = "DAPI"
UNMIXED_DAPI_KEY = "dapi_unmixed"
CELLPOSE_KEY = "cellpose_boundaries"
STARDIST_KEY = "stardist_boundaries"


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Run nuclei segmentation in Xenium with optional DAPI unmixing and background removal"
    )

    def _existing_path(p: str) -> Path:
        pth = Path(p)
        if not pth.exists():
            raise argparse.ArgumentTypeError(f"Path does not exist: {pth}")
        return pth

    parser.add_argument(
        "--xenium-file",
        dest="xenium_file",
        required=True,
        type=_existing_path,
        help="Input SpatialData Zarr produced from Xenium.",
    )
    parser.add_argument(
        "--out", required=True, type=Path, help="Output Zarr path to write (e.g., <sample>_reseg_cells.zarr)."
    )
    parser.add_argument(
        "--segment-config",
        dest="segment_config",
        required=True,
        type=_existing_path,
        help="JSON file with method-specific segmentation parameters.",
    )
    parser.add_argument(
        "--patch-config",
        dest="patch_config",
        required=True,
        type=_existing_path,
        help="JSON file with image patching parameters.",
    )

    parser.add_argument("--n-workers", default=1, type=int, required=False, help="Number of Dask workers (default: 1).")
    parser.add_argument(
        "--unmix",
        action="store_true",
        help=f"Compute an unmixed '{DAPI_KEY}' channel as '{UNMIXED_DAPI_KEY}' and use it.",
    )
    parser.add_argument(
        "--method",
        required=True,
        choices=["cellpose", "stardist"],
        help="Segmentation method to use.",
    )

    return parser.parse_args()


def setup_sopa(n_workers: int):
    sopa.settings.parallelization_backend = "dask"
    sopa.settings.dask_client_kwargs = {
        "n_workers": n_workers,
        "processes": False,
        "threads_per_worker": 1,
        "dashboard_address": ":8787",
    }
    sopa.log.setLevel(sopa.logging.ERROR)


def subtract_background(img, radius, method="rolling_ball"):
    if method == "rolling_ball":
        bg = rolling_ball(img, radius=radius)
        fg = np.maximum(img - bg, 0)
    elif method == "white_tophat":
        footprint = disk(int(max(1, radius)))
        fg = white_tophat(img, footprint=footprint)
        bg = img - fg
    else:
        raise ValueError(f"Unknown bg method: {method!r}")
    return fg, bg


def stratified_bg_sample(
    y: np.ndarray,
    bg_nonzero: np.ndarray,
    max_samples: int = 2_000_000,
    n_bins: int = 32,
    seed: int = 0,
) -> np.ndarray:
    if bg_nonzero.size <= max_samples:
        return bg_nonzero

    rng = np.random.default_rng(seed)
    bg_values = y[bg_nonzero]

    edges = np.quantile(bg_values, np.linspace(0, 1, n_bins + 1))
    edges = np.unique(edges)

    if edges.size < 3:
        return rng.choice(bg_nonzero, size=max_samples, replace=False)

    bin_id = np.digitize(bg_values, edges[1:-1], right=True)
    per_bin = max_samples // (edges.size - 1)

    chosen = []
    for b in range(edges.size - 1):
        candidates = bg_nonzero[bin_id == b]
        if candidates.size == 0:
            continue
        take = min(per_bin, candidates.size)
        chosen.append(rng.choice(candidates, size=take, replace=False))

    if not chosen:
        return rng.choice(bg_nonzero, size=max_samples, replace=False)

    sampled = np.concatenate(chosen)
    return sampled


def regress_out_channel(
    target_channel,
    spillover_channel,
    keep_baseline: bool = False,
    bg_percentile: float = 35.0,
    rescale: bool = True,
):
    h, w = target_channel.shape
    features = spillover_channel.reshape(1, -1).T

    y = target_channel.reshape(-1)
    X = np.concatenate(
        [np.ones((h * w, 1)), features],
        axis=1,
    )

    bg_threshold = np.percentile(target_channel, bg_percentile)
    bg_mask = y <= bg_threshold
    bg_nonzero = np.flatnonzero(bg_mask)

    if bg_nonzero.size == 0:
        raise ValueError("No background pixels selected; try increasing bg_percentile.")

    bg_idx = stratified_bg_sample(y, bg_nonzero, n_bins=8)

    beta, _ = nnls(X[bg_idx], y[bg_idx])
    if keep_baseline:
        weights = beta[1:, 0]
        predicted_spillover = (features @ weights).reshape(h, w)
    else:
        predicted_spillover = (X @ beta).reshape(h, w)

    predicted_spillover = np.minimum(predicted_spillover, target_channel)
    target_unmixed_f64 = target_channel - predicted_spillover

    p_min = 0
    p_max = np.percentile(target_unmixed_f64, 99.5)

    if rescale:
        u16_max = np.iinfo(np.uint16).max
        target_unmixed = rescale_intensity(
            target_unmixed_f64,
            in_range=(p_min, p_max),
            out_range=(0, u16_max),
        ).astype(np.uint16)
    else:
        target_unmixed = np.clip(target_unmixed_f64, 0, None).astype(np.uint16)

    return target_unmixed, predicted_spillover, beta


def unmix_channel(
    image,
    target_channel: str = "DAPI",
    spillover_channel: str = "ATP1A1/CD45/E-Cadherin",
    keep_baseline: bool = False,
    bg_percentile: float = 35.0,
    bg_method: str = "rolling_ball",
    rb_radius: int = 9,
    subtract_bg: bool = False,
    rescale: bool = True,
):
    """
    Unmix `target_channel` by regressing out spillover from `spillover_channel`.

    Parameters
    ----------
    image : SpatialData Image2DModel object
    target_channel : str
    spillover_channel : str
    keep_baseline : bool
        Passed to regress_out_channel.
    bg_percentile : float
        Passed to regress_out_channel.
    bg_method : str
        Used only if subtract_bg=True.
    rb_radius : int
        Used only if subtract_bg=True.
    subtract_bg : bool
        If True, run subtract_background(...) on both target and spillover before regression.

    Returns
    -------
    unmixed_element : spatialdata Image2DModel
    """

    image_data = image["scale0"]["image"]
    transformation = spatialdata.transformations.get_transformation(image)

    target_u16_raw = image_data.sel(c=target_channel).values
    spillover_u16_raw = image_data.sel(c=spillover_channel).values

    if subtract_bg:
        target_u16, _ = subtract_background(target_u16_raw, rb_radius, bg_method)
        spillover_u16, _ = subtract_background(spillover_u16_raw, rb_radius, bg_method)
    else:
        target_u16 = target_u16_raw
        spillover_u16 = spillover_u16_raw

    target_unmixed_u16, predicted_spillover, beta = regress_out_channel(
        target_channel=target_u16,
        spillover_channel=spillover_u16,
        keep_baseline=keep_baseline,
        bg_percentile=bg_percentile,
        rescale=rescale,
    )

    arr_cyx = target_unmixed_u16[None, ...]
    unmixed_element = spatialdata.models.Image2DModel.parse(
        data=arr_cyx,
        dims=("c", "y", "x"),
        c_coords=[target_channel],
        transformations={"global": transformation},
    )

    return unmixed_element


def main():
    args = _parse_args()
    setup_sopa(args.n_workers)

    sdata = spatialdata.read_zarr(args.xenium_file.resolve())
    segment_config = json.loads(args.segment_config.read_bytes())
    patches_config = json.loads(args.patch_config.read_bytes())

    if args.unmix:
        image_key = sdata.attrs["xenflow"]["current"]["morphology_image"]
        dapi_unmixed = unmix_channel(sdata.images[image_key], DAPI_KEY)
        sdata.images[UNMIXED_DAPI_KEY] = dapi_unmixed
        sdata.attrs["xenflow"]["current"]["morphology_image"] = UNMIXED_DAPI_KEY

    image_key = sdata.attrs["xenflow"]["current"]["morphology_image"]
    warnings.warn(f"{image_key=}")

    sopa.make_image_patches(sdata, image_key=image_key, **patches_config)

    if args.method == "cellpose":
        sopa.segmentation.cellpose(
            sdata,
            image_key=image_key,
            channels=[DAPI_KEY],
            key_added=CELLPOSE_KEY,
            **segment_config,
        )
        sdata.attrs["xenflow"]["current"]["nucleus_shapes"] = CELLPOSE_KEY
    elif args.method == "stardist":
        sopa.segmentation.stardist(
            sdata,
            image_key=image_key,
            channels=[DAPI_KEY],
            key_added=STARDIST_KEY,
            **segment_config,
        )
        sdata.attrs["xenflow"]["current"]["nucleus_shapes"] = STARDIST_KEY

    sdata.write(args.out.resolve())


if __name__ == "__main__":
    main()
