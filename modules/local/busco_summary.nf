process BUSCO_SUMMARY {

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path(summaries_domain)
    path(summaries_specific)
    path(failed_bins)

    output:
    path ("*-busco_summary.tsv"), emit: summary

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def auto = params.busco_reference ? "" : "-a"
    def ss = summaries_specific.sort().size() > 0 ? "-ss ${summaries_specific}" : ""
    def sd = summaries_domain.sort().size() > 0 ? "-sd ${summaries_domain}" : ""
    def f = ""
    if (!params.busco_reference && failed_bins.sort().size() > 0)
        f = "-f ${failed_bins}"
    """
    summary_busco.py $auto $ss $sd $f -o "${task.ext.prefix}-busco_summary.tsv"
    """
}