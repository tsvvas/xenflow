/*
 * Module: DETECT_TISSUE
 * Detects tissue contours and bounding boxes in Xenium data
 */
process DETECT_TISSUE {
    tag        "${sample_id}"
    publishDir "${workflow.launchDir}/${params.outdir_root}/${params.outdir_coords}", mode: 'copy'
    container  "${params.containerdir}/sopa.sif"

    cpus    params.cpus
    memory  params.mem
    time    params.time
    
    beforeScript '''
    export NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
    export MPLCONFIGDIR=${XDG_CACHE_HOME}/.matplotlib
    '''

    input:
        tuple path(zarr_file), val(sample_id)

    output:
        tuple \
            path("${sample_id}_bbox_coords.json"), \
            path("${sample_id}_regions.zarr"), \
            path("${sample_id}_coords.png"), \
            val(sample_id)

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spatial

    get_tissue_boundaries.py \
        --xenium-file  ${zarr_file} \
        --kernel-size  ${params.detect_tissue_kernel_size} \
        --out-bbox     ${sample_id}_bbox_coords.json \
        --out-regions  ${sample_id}_regions.zarr \
        --out-fig      ${sample_id}_coords.png \
    """
}
