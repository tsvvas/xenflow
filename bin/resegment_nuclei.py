#!/usr/bin/env python
import argparse
import json
from pathlib import Path

import numpy as np
import sopa
import spatialdata
from scipy.optimize import nnls
from skimage.exposure import rescale_intensity
from skimage.morphology import disk, opening, white_tophat
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
        bg = opening(img, footprint=footprint)
        fg = white_tophat(img, footprint=footprint)
    else:
        raise ValueError(f"Unknown bg method: {method!r}")
    return fg, bg


def unmix_channel(
    image,
    target_channel: str = "DAPI",
    spillover_channels: list[str] | None = None,
    keep_baseline: bool = False,
    bg_percentile: float = 35.0,
    bg_method="rolling_ball",
    rb_radius: int = 120,
):
    if spillover_channels is None:
        spillover_channels = ["ATP1A1/CD45/E-Cadherin"]

    image_data = image["scale0"]["image"]
    transformation = spatialdata.transformations.get_transformation(image)

    target_f64 = image_data.sel(c=target_channel).astype(np.float64).values
    target_f64, _ = subtract_background(target_f64, rb_radius, bg_method)
    spillover_f64_raw = [image_data.sel(c=ch).astype(np.float64).values for ch in spillover_channels]
    spillover_f64 = np.stack([subtract_background(ch, rb_radius, bg_method)[0] for ch in spillover_f64_raw], axis=0)

    h, w = target_f64.shape
    k = spillover_f64.shape[0]
    features = spillover_f64.reshape(k, -1).T

    y = target_f64.reshape(-1)
    X = np.concatenate(
        [np.ones((h * w, 1), dtype=np.float64), features],
        axis=1,
    )

    thr = np.percentile(target_f64, bg_percentile)
    idx = (target_f64 <= thr).ravel()
    beta, _ = nnls(X[idx], y[idx])
    if keep_baseline:
        weights = beta[1:, 0]
        predicted_spillover = (features @ weights).reshape(h, w)
    else:
        predicted_spillover = (X @ beta).reshape(h, w)

    predicted_spillover = np.minimum(predicted_spillover, target_f64)
    target_unmixed_f64 = target_f64 - predicted_spillover

    positive_values = target_unmixed_f64[target_unmixed_f64 > 0]
    if positive_values.size >= 16:
        p_min, p_max = np.percentile(positive_values, (0.5, 99.5))
    else:
        p_min, p_max = np.percentile(target_unmixed_f64, (0.5, 99.5))
    u16_max = np.iinfo(np.uint16).max
    epsilon = 1
    target_unmixed_u16 = rescale_intensity(
        target_unmixed_f64,
        in_range=(p_min, p_max),
        out_range=(epsilon, u16_max),
    ).astype(np.uint16)

    arr_cyx = target_unmixed_u16[None, ...]

    unmixed_element = spatialdata.models.Image2DModel.parse(
        data=arr_cyx,
        dims=("c", "y", "x"),
        c_coords=[target_channel],
        transformations={"global": transformation},
        scale_factors=None,
        chunks=(1, min(512, h), min(512, w)),
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
    else:
        image_key = sdata.attrs["xenflow"]["current"]["morphology_image"]

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
