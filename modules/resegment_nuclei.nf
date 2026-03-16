/*
 * Module: RESEGMENT_NUCLEI
 * Re-segments nuclei from Xenium images (Cellpose or StarDist) and writes
 * updated labels into a SpatialData-formatted Zarr.
 */

process RESEGMENT_NUCLEI {
    tag "${sample_id}"
    publishDir "${workflow.launchDir}/${params.outdir_root}/${params.outdir_reseg_nuclei}", mode: 'copy'
    container "${params.containerdir}/sopa_gpu.sif"

    cpus params.cpus
    memory params.mem
    time params.time
    queue params.queue
    clusterOptions params.cluster_opts

    beforeScript '''
    export NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
    export MPLCONFIGDIR=${XDG_CACHE_HOME}/.matplotlib
    '''

    input:
    tuple path(zarr_file), val(sample_id)

    output:
    tuple path("${sample_id}_reseg_nuclei.zarr"), val(sample_id)

    script:
    def method = params.nuclei_segment_method?.toLowerCase()
    if (method != 'cellpose' && method != 'stardist') {
        throw new IllegalArgumentException(
            "Unknown nuclei segmentation method: '${params.nuclei_segment_method}'. Supported: 'cellpose', 'stardist'."
        )
    }

    def segment_config_json = (method == 'cellpose'
        ? JsonUtils.toPrettyJson(params.cellpose)
        : JsonUtils.toPrettyJson(params.stardist))
    def patch_config_json = JsonUtils.toPrettyJson(params.image_patches)
    def unmix_flag = params.unmix_dapi ? '--unmix' : ''

    """
    export CELLPOSE_LOCAL_MODELS_PATH=\$XDG_CACHE_HOME/.cellpose
    export DASK_DATAFRAME__QUERY_PLANNING=False
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spatial

    printf '%s\n' '${segment_config_json}' > segment_config.json
    printf '%s\n' '${patch_config_json}' > patch_config.json

    resegment_nuclei.py \
        --xenium-file     ${zarr_file} \
        --out             ${sample_id}_reseg_nuclei.zarr \
        --method          ${method} \
        --segment-config  segment_config.json \
        --patch-config    patch_config.json \
        --n-workers       1 \
        ${unmix_flag}
    """
}
