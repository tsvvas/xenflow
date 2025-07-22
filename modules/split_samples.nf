/*
 * Module: SPLIT_SAMPLES
 * Adds region polygons to the dataset, splits the SpatialData object
 * and writes one AnnData table (`.h5ad`) per region.
 */
process SPLIT_SAMPLES {

    tag        "${sample_id}"
    publishDir "${workflow.launchDir}/${params.outdir_h5ad}", mode: 'copy'
    container  "${params.containerdir}/sopa.sif"   

    cpus           params.cpus
    memory         params.mem
    time           params.time
    queue          params.queue
    clusterOptions params.cluster_opts

    input:
        tuple path(zarr_file), path(regions_zarr), val(sample_id)

    output:
        tuple path("${sample_id}_*.h5ad"), val(sample_id)

    beforeScript '''
    export NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
    export MPLCONFIGDIR=${XDG_CACHE_HOME}/.matplotlib
    '''

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spatial

    split_samples_by_region.py \
        --dataset-zarr  ${zarr_file} \
        --regions-zarr  ${regions_zarr} \
        --sample-id     ${sample_id} \
    """
}
