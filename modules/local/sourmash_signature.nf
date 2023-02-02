/*
-----------------------------Sourmash added for Taxonomic Classification by Zifo----------------------------
*/
process SOURMASH_SIGNATURE {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::sourmash=4.4.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.4.0--hdfd78af_0' :
        'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path('*.sig')	    , emit: signatures
    path "versions.yml"                 , emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
    """
    sourmash sketch dna ${args} ${bins} --output "${prefix}.sig"   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version)
    END_VERSIONS
    """
}