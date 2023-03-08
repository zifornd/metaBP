process MAG_DEPTHS_SUMMARY {

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path(mag_depths)

    output:
    path("bin_depths_summary-*.tsv"), emit: summary

    script:
    def args = task.ext.args ?: ''
    """
    get_mag_depths_summary.py --depths ${mag_depths} \
                              --out "bin_depths_summary-${task.ext.prefix}.tsv"
    """
}