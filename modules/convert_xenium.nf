/*
 * Module: CONVERT_XENIUM
 * Converts a Xenium run to a Zarr store with SpatialData-IO
 */
process CONVERT_XENIUM {                         
    tag          "${sample_id}"
    publishDir   "${workflow.launchDir}/${params.outdir_root}/${params.outdir_zarr}", mode: 'copy'
    container    "${params.containerdir}/sopa.sif"

    cpus           params.cpus
    memory         params.mem
    time           params.time
    queue          params.queue
    clusterOptions params.cluster_opts

    input:
        tuple path(xenium_dir), val(sample_id)

    output:
        tuple path("${sample_id}.zarr"), val(sample_id)

    beforeScript '''
    export NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
    export MPLCONFIGDIR=${XDG_CACHE_HOME}/.matplotlib
    '''

    script:
    """
    export RAPIDS_NO_INITIALIZE=1
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spatial

    xenium_to_zarr.py \
        --xenium_dir ${xenium_dir} \
        --out        ${sample_id}.zarr
    """
}
