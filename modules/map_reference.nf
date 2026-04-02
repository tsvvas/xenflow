/*
 * Module: MAP_REFERENCE
 * Maps cell types from reference single-cell dataset to the spatial.
 */

process MAP_REFERENCE {
    tag        "${sample_id}"
    publishDir "${workflow.launchDir}/${params.outdir_root}/${params.outdir_mapped}", mode: 'copy'
    container  "${params.containerdir}/sopa_gpu.sif" 

    cpus           params.cpus
    memory         params.mem
    time           params.time
    queue          params.queue
    clusterOptions params.cluster_opts

    input:
        tuple path(h5ad_file), val(sample_id)

    output:
        tuple \
            path("${h5ad_file.baseName}_mapped.h5ad"), \
            path("${sample_id}_training_scores.png"), \
            val(sample_id)

    script:
    """
    export uid=${h5ad_file.baseName}
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spatial

    map_reference.py \
        --xenium-file ${h5ad_file} \
        --reference-file ${params.reference_dataset} \
        --out \${uid}_mapped.h5ad \
        --cell-type-key ${params.cell_type_key} \
        --training-plot ${sample_id}_training_scores.png
    """

}
