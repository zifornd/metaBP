process MAG_DEPTHS {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"

    // Using container from metabat2 process, since this will be anyway already downloaded and contains biopython and pandas
    conda (params.enable_conda ? "bioconda::metabat2=2.15 conda-forge::python=3.6.7 conda-forge::biopython=1.74 conda-forge::pandas=1.1.5" : null)
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' :
        'quay.io/biocontainers/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' }"

    input:
    tuple val(meta), path(bins)
    path(contig_depths)

    output:
    tuple val(meta), path("${meta.assembler}-${meta.trimmer}-${meta.id}-binDepths.tsv"), emit: depths

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
    """
    get_mag_depths.py --bins ${bins} \
                    --depths ${contig_depths} \
                    --assembly_name "${prefix}" \
                    --out "${prefix}-binDepths.tsv"
    """
}