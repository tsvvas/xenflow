/*
 * Module: RESEGMENT_CELLS
 * Re-segments cells using transcript molecules (e.g., Baysor or ProSeg)
 * and writes updated labels into a SpatialData-formatted Zarr.
 */
import static groovy.json.JsonOutput.prettyPrint
import static groovy.json.JsonOutput.toJson

process RESEGMENT_CELLS {                         
    tag          "${sample_id}"
    publishDir   "${workflow.launchDir}/${params.outdir_root}/${params.outdir_reseg}", mode: 'copy'
    container    "${params.containerdir}/sopa.sif"

    cpus           params.cpus
    memory         params.mem
    time           params.time
    queue          params.queue
    clusterOptions params.cluster_opts

    input:
        tuple path(zarr_file), val(sample_id)

    output:
        tuple path("${sample_id}_reseg_cells.zarr"), val(sample_id)

    beforeScript '''
    export NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
    export MPLCONFIGDIR=${XDG_CACHE_HOME}/.matplotlib
    '''

    script:
    def segment_config_json
    switch (params.cell_segment_method?.toLowerCase()) {
        case 'baysor':
            segment_config_json = prettyPrint(toJson(params.baysor))
            break
        case 'proseg':
            // Proseg wrapper in sopa doesn't have config
            segment_config_json = ""
            break
        default:
            throw new IllegalArgumentException(
                "Unknown cell segmentation method: '${params.cell_segment_method}'. " +
                "Supported: 'baysor', 'proseg'."
            )
    }

    def patch_config_json = prettyPrint(toJson(params.transcript_patches))

    """
    export RAPIDS_NO_INITIALIZE=1
    export DASK_DATAFRAME__QUERY_PLANNING=False
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spatial

    printf '%s\n' '${segment_config_json}' > segment_config.json
    printf '%s\n' '${patch_config_json}' > patch_config.json

    resegment_cells.py \
        --xenium_file     ${zarr_file} \
        --out             ${sample_id}_reseg_cells.zarr \
        --method          ${params.cell_segment_method} \
        --segment-config  segment_config.json \
        --patch-config    patch_config.json \
        --n-workers       ${task.cpus}
    """
}
