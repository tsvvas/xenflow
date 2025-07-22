#!/usr/bin/env python
import argparse
import json
from pathlib import Path
from typing import List, Optional, Tuple

import cv2
import geopandas as gpd
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import spatialdata
from shapely.geometry import Polygon
from spatialdata.models import ShapesModel


class XeniumTissueDetector:
    """
    A class to detect contours and their respective coordinates from Xenium slide images using cell labels.

    Attributes:
    ----------
    bimage : np.ndarray
        Binary image used for contour detection.
    aff_mtx : np.ndarray
        Affine transformation matrix to convert coordinates from pixel coordinate system to physical space.
    contours : Optional[List[np.ndarray]]
        List of contours detected in the binary image.
    transformed_contours : List[np.ndarray]
        List of transformed contours, with coordinates in physical space.
    sample_coordinates : List[Tuple[float, float, float, float]]
        List of transformed bounding box coordinates, stored as (x1, y1, x2, y2).
    """

    def __init__(self, sdata) -> None:
        """
        Initializes the XeniumTissueDetector object with transformations based on spatial data.

        Parameters:
        ----------
        sdata : object
            Spatial data containing cell labels and shape transformations for computing the affine matrix.
        """
        self.bimage: np.ndarray = (
            sdata.labels["cell_labels"]["scale4"].image.data.clip(min=0, max=1).compute().astype(np.uint8)
        )

        mat1 = spatialdata.transformations.get_transformation(sdata.shapes["cell_boundaries"]).to_affine_matrix(
            input_axes=("x", "y"), output_axes=("x", "y")
        )
        mat2 = spatialdata.transformations.get_transformation(
            sdata.labels["cell_labels"]["scale4"].image
        ).to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))

        self.aff_mtx: np.ndarray = np.linalg.inv((np.linalg.inv(mat2) @ mat1))
        self.tf2global = spatialdata.transformations.get_transformation(sdata.shapes["cell_boundaries"])
        self.contours: Optional[List[np.ndarray]] = None
        self.transformed_contours: List[np.ndarray] = []
        self.sample_coordinates: List[Tuple[float, float, float, float]] = []
        self.sdata_regions: spatialdata.SpatialData | None = None

    def apply_affine_transformation(self, points: np.ndarray) -> np.ndarray:
        """
        Applies the affine transformation matrix to a set of points.

        Parameters:
        ----------
        points : np.ndarray
            Array of points to be transformed, with shape (n_points, 2).

        Returns:
        -------
        np.ndarray
            Transformed points with shape (n_points, 2).
        """
        points_homogeneous = np.hstack((points, np.ones((points.shape[0], 1))))
        return (points_homogeneous @ self.aff_mtx.T)[:, :2]

    def find_contours(self, kernel_size: int = 150) -> None:
        """
        Detects contours in the binary image using dilation and morphological closing operations.

        This method applies transformations to detected contours and their corresponding bounding boxes.

        Parameters:
        ----------
        kernel_size : int, optional
            Size of the morphological kernel used for contour merging. Default is 150.
        """
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size, kernel_size))
        dilated = cv2.dilate(self.bimage, kernel, iterations=2)
        closed = cv2.morphologyEx(dilated, cv2.MORPH_CLOSE, kernel)

        # Find contours in the processed image
        self.contours, _ = cv2.findContours(closed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        # Process each contour
        for contour in self.contours:
            # Transform and store contour points
            transformed_contour = self.apply_affine_transformation(contour.reshape(-1, 2))
            self.transformed_contours.append(transformed_contour)

            # Calculate and transform bounding box
            x, y, w, h = cv2.boundingRect(contour)
            x1, y1, x2, y2 = x, y, x + w, y + h
            x1, y1, x2, y2 = transform_bounding_rects([[x1, y1, x2, y2]], self.aff_mtx)[0]
            self.sample_coordinates.append((x1, y1, x2, y2))

    def collect_regions_to_spatialdata(self):
        gdf = gpd.GeoDataFrame(
            {
                "region_id": [f"region{i + 1:02d}" for i, _ in enumerate(self.transformed_contours)],
                "geometry": [Polygon(pix_xy) for pix_xy in self.transformed_contours],
            }
        )

        regions_shapes = ShapesModel.parse(gdf, transformations={"global": self.tf2global})
        self.sdata_regions = spatialdata.SpatialData()
        self.sdata_regions.shapes["regions"] = regions_shapes

    def plot_contours(
        self,
        contour_color: str = "blue",
        rect_color: str = "black",
        rect_linewidth: int = 2,
        **kwargs,
    ) -> plt.Axes:
        """
        Plots the transformed contours and their bounding boxes on a 2D plot.

        Parameters:
        ----------
        contour_color : str, optional
            Color used to plot the contour points. Default is "blue".
        rect_color : str, optional
            Color used to draw the bounding rectangles. Default is "black".
        rect_linewidth : int, optional
            Line width of the bounding rectangles. Default is 2.
        **kwargs : dict, optional
            Additional keyword arguments passed to `matplotlib.pyplot.subplots`.

        Returns:
        -------
        matplotlib.axes.Axes
            Axes object containing the plot.

        Raises:
        ------
        ValueError
            If `find_contours` has not been called and no contours have been detected.
        """
        # Check if contours have been found and transformed
        if not self.transformed_contours or not self.sample_coordinates:
            raise ValueError("Contours not found. Please run `find_contours` first.")

        fig, ax = plt.subplots(**kwargs)

        # Plot contours and rectangles
        for transformed_contour, (x1, y1, x2, y2) in zip(self.transformed_contours, self.sample_coordinates):
            ax.scatter(transformed_contour[:, 0], transformed_contour[:, 1], c=contour_color)
            ax.add_patch(
                matplotlib.patches.Rectangle(
                    (x1, y1),
                    x2 - x1,
                    y2 - y1,
                    edgecolor=rect_color,
                    fill=False,
                    lw=rect_linewidth,
                )
            )
        return ax


def transform_bounding_rects(bounding_rects: list, affine_matrix: np.ndarray) -> list:
    """
    Transforms a list of bounding rectangles using the provided affine transformation matrix.

    Args:
        bounding_rects (list of tuples): List of bounding rectangles represented as [(x1, y1, x2, y2), ...], where:
            - x1 (float): Top-left x-coordinate of the rectangle.
            - y1 (float): Top-left y-coordinate of the rectangle.
            - x2 (float): Bottom-right x-coordinate of the rectangle.
            - y2 (float): Bottom-right y-coordinate of the rectangle.
        affine_matrix (np.ndarray): 3x3 affine transformation matrix used to transform the rectangles.

    Returns:
        list of tuples: Transformed bounding rectangles represented as [(new_x1, new_y1, new_x2, new_y2), ...], where:
            - new_x1 (float): New top-left x-coordinate of the transformed rectangle.
            - new_y1 (float): New top-left y-coordinate of the transformed rectangle.
            - new_x2 (float): New bottom-right x-coordinate of the transformed rectangle.
            - new_y2 (float): New bottom-right y-coordinate of the transformed rectangle.
    """
    transformed_rects = []
    for x1, y1, x2, y2 in bounding_rects:
        corners = np.array([[x1, y1, 1], [x2, y1, 1], [x1, y2, 1], [x2, y2, 1]])
        transformed_corners = (affine_matrix @ corners.T).T

        transformed_x = transformed_corners[:, 0]
        transformed_y = transformed_corners[:, 1]

        new_x1 = transformed_x.min()
        new_y1 = transformed_y.min()
        new_x2 = transformed_x.max()
        new_y2 = transformed_y.max()

        transformed_rects.append((new_x1, new_y1, new_x2, new_y2))

    return transformed_rects


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--xenium-file", required=True, type=Path)
    p.add_argument("--kernel-size", required=False, type=int, default=150)
    p.add_argument(
        "--out-bbox",
        required=True,
        type=Path,
        help="where to write the sample-coordinates JSON",
    )
    p.add_argument("--out-regions", required=True, type=Path)
    p.add_argument("--out-fig", required=True, type=Path, help="PNG image of the contours")
    args = p.parse_args()

    sdata = spatialdata.read_zarr(args.xenium_file.resolve())

    xeniumcv = XeniumTissueDetector(sdata)
    xeniumcv.find_contours(kernel_size=args.kernel_size)
    xeniumcv.collect_regions_to_spatialdata()

    xeniumcv.sdata_regions.write(args.out_regions.resolve())

    ax = xeniumcv.plot_contours(figsize=(5, 8))
    ax.xaxis.tick_top()
    plt.gca().invert_yaxis()
    plt.savefig(args.out_fig, dpi=600, bbox_inches="tight")

    sample_coordinates = {
        "sample_coords": {
            f"region{i + 1:0>2}": xeniumcv.sample_coordinates[i] for i in range(len(xeniumcv.sample_coordinates))
        }
    }

    with Path(args.out_bbox.resolve()).open("w", encoding="UTF-8") as target:
        json.dump(sample_coordinates, target, indent=4)


if __name__ == "__main__":
    main()
