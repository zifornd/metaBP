
process MAG_DEPTHS_PLOT {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas=1.3.0 anaconda::seaborn=0.11.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0' :
        'quay.io/biocontainers/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0' }"

    input:
    tuple val(meta), path(depths)
    path(sample_groups)

    output:
    tuple val(meta), path("${meta.assembler}-${meta.trimmer}-${meta.id}-binDepths.heatmap.png"), emit: heatmap

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
    """
    plot_mag_depths.py --bin_depths ${depths} \
                    --groups ${sample_groups} \
                    --out "${prefix}-binDepths.heatmap.png"
    """
}
