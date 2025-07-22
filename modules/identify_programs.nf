/*
 * Module: IDENTIFY_PROGRAMS
 * Runs consensus-NMF on one ROI and keeps **all** output.
 */
process IDENTIFY_PROGRAMS {

    tag        "${h5ad_file.baseName}"                     
    publishDir "${workflow.launchDir}/${params.outdir_programs}", mode: 'copy'
    container  "${params.containerdir}/sopa.sif"

    cpus    params.cpus
    memory  params.mem
    time    params.time

    input:
        path h5ad_file

    output:
        path "${h5ad_file.baseName}_cnmf"

    script:
    """
    export uid=${h5ad_file.baseName}
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spatial

    identify_cnmf_programs.py \
        --h5ad              ${h5ad_file} \
        --out-dir           \${uid}_cnmf \
        --components        ${params.cnmf_components} \
        --n-iter            ${params.cnmf_n_iter} \
        --num-highvar-genes ${params.cnmf_n_highvar} \
        --seed              ${params.seed} \
        --workers           ${task.cpus}
    """
}
