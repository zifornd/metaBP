process QUAST_SPADES {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"

    conda (params.enable_conda ? "bioconda::quast=5.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2' :
        'quay.io/biocontainers/quast:5.0.2--py37pl526hb5aa323_2' }"

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(contig)

    output:
    path "QUAST/*"          , emit: qc
    path 'versions.yml'     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
       """
       metaquast.py --threads "${task.cpus}" --max-ref-number 0 -l "${prefix}_scaffolds","${prefix}_contigs" "${assembly}" "${contig}" -o "QUAST"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        metaquast: \$(metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//")
    END_VERSIONS
    """
}