/*
--------------------------Sourmash added for Taxonomic Classification by Zifo----------------------
*/
process SOURMASH_SUMMARIZE {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"

    conda (params.enable_conda ? 'bioconda::sourmash=4.4.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.4.0--hdfd78af_0' :
        'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0' }"

    input:
    path(db), stageAs: 'database.gz'
    tuple val(meta), path(signatures)

    output:
    tuple val(meta), path('*.csv')	    , emit: report
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def ram = task.memory
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
    """
    echo ${ram} > mem.txt
    sourmash lca summarize \\
        --db ${db} \\
        --query ${signatures} \\
        -o '${prefix}.csv' \
        2> sourmash.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version)
    END_VERSIONS
    """
}