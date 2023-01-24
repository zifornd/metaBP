process QUAST_BINS_SUMMARY {

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"
    
    input:
    path(summaries)

    output:
    path("*-quast_summary.tsv"), emit: summary

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    QUAST_BIN=\$(echo \"$summaries\" | sed 's/[][]//g')
    IFS=', ' read -r -a quast_bin <<< \"\$QUAST_BIN\"
    for quast_file in \"\${quast_bin[@]}\"; do
        if ! [ -f "${prefix}-quast_summary.tsv" ]; then
            cp "\${quast_file}" "${prefix}-quast_summary.tsv"
        else
            tail -n +2 "\${quast_file}" >> "${prefix}-quast_summary.tsv"
        fi
    done
    """
}
