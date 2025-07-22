# Xenflow

Xenium spatial transcriptomics analysis pipeline built with Nextflow and Singularity/Apptainer. The pipeline uses cv2 to automatically identify tissue contours and bounding boxes; it is designed to analyse tissue microarrays and multi-sample slides.

## Quick start
1. Build a container `sopa.sif` from the definition file from [this repository](https://github.com/tsvvas/singularity_defs).
2. Export `CONTAINERDIR` and `PROJECTDIR` variables. The former is used to find the container, and the latter to mount the data to singularity during run.
3. Adjust resources in [resources.config](config/resources.config) for pipeline steps.
4. Create `run01.config` in config directory from [run.template](config/run.template) to configure step parameters.
5. Create environment with nextflow installation using `make env`.
6. Run the pipeline using `make run`.

## Key steps
The pipeline has 4 steps starting from raw xenium output to identification of gene programs using [cNMF](https://github.com/dylkot/cNMF).

- CONVERT_XENIUM converts raw data to [spatialdata](https://github.com/scverse/spatialdata)-formatted zarr archive for downstream processing.
- DETECT_TISSUE automatically identifies tissue contours and bouding boxes to split multi-sample slides into separate sample objects for independent analysis. Especially helpful when working with tissue microarrays.
- SPLIT_SAMPLES creates one AnnData h5ad archive per sample based on the tissue contours.
- IDENTIFY_PROGRAMS runs cNMF for each sample. Note that it doesn't automatically select the best number of programs in this pipeline leaving this task for user.
