process BIN_SUMMARY {

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path(bin_depths)
    path(busco_sum)
    path(quast_sum)
    path(gtdbtk_sum)

    output:
    path("bin_summary_*.tsv"), emit: summary

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def busco_summary  = busco_sum.sort().size() > 0 ?  "--busco_summary ${busco_sum}" : ""
    def quast_summary  = quast_sum.sort().size() > 0 ?  "--quast_summary ${quast_sum}" : ""
    def gtdbtk_summary = gtdbtk_sum.sort().size() > 0 ? "--gtdbtk_summary ${gtdbtk_sum}" : ""
    """
    combine_tables.py --depths_summary ${bin_depths} \
                      $busco_summary \
                      $quast_summary \
                      $gtdbtk_summary \
                      --out bin_summary_${prefix}.tsv
    """
}