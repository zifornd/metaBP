process GTDBTK_SUMMARY {

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path(qc_discarded_bins)
    path(gtdbtk_summaries)
    path(filtered_bins)
    path(failed_bins)

    output:
    path ("gtdbtk_summary-*.tsv"), emit: summary

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def discarded = qc_discarded_bins.sort().size() > 0 ? "--qc_discarded_bins ${qc_discarded_bins}" : ""
    def summaries = gtdbtk_summaries.sort().size() > 0 ?  "--summaries ${gtdbtk_summaries}" : ""
    def filtered  = filtered_bins.sort().size() > 0 ?     "--filtered_bins ${filtered_bins}" : ""
    def failed    = failed_bins.sort().size() > 0 ?       "--failed_bins ${failed_bins}" : ""
    """
    summary_gtdbtk.py $args $discarded $summaries $filtered $failed --out ${prefix}-gtdbtk_summary.tsv
    """
}