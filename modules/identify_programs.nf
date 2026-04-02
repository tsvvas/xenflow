/*
 * Module: IDENTIFY_PROGRAMS
 * Runs consensus-NMF on one ROI and keeps **all** output.
 */
process IDENTIFY_PROGRAMS {

    tag        "${h5ad_file.baseName}"                     
    publishDir "${workflow.launchDir}/${params.outdir_root}/${params.outdir_programs}", mode: 'copy'
    container  "${params.containerdir}/sopa.sif"

    cpus    params.cpus
    memory  params.mem
    time    params.time

    input:
        tuple path(h5ad_file), val(sample_id)

    output:
        path "${sample_id}_cnmf"

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spatial

    identify_cnmf_programs.py \
        --h5ad              ${h5ad_file} \
        --id                ${sample_id} \
        --out-dir           ${sample_id}_cnmf \
        --components        ${params.cnmf_components} \
        --n-iter            ${params.cnmf_n_iter} \
        --num-highvar-genes ${params.cnmf_n_highvar} \
        --seed              ${params.seed} \
        --workers           ${task.cpus}
    """
}
