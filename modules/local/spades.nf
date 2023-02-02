process SPADES {
    tag "${meta.trimmer}-${meta.id}"

    conda (params.enable_conda ? "bioconda::spades=3.15.3" : null)
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.3--h95f258a_0' :
        'quay.io/biocontainers/spades:3.15.3--h95f258a_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.trimmer}-${meta.id}_scaffolds.fasta"), emit: assembly
    tuple val(meta), path("${meta.trimmer}-${meta.id}_contigs.fasta")  , emit: contig
    path "${meta.trimmer}-${meta.id}.log"                              , emit: log
    path "${meta.trimmer}-${meta.id}_contigs.fasta.gz"                 , emit: contigs_gz
    path "${meta.trimmer}-${meta.id}_scaffolds.fasta.gz"               , emit: assembly_gz
    path "${meta.trimmer}-${meta.id}_graph.gfa.gz"                     , emit: graph
    path 'versions.yml'                                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix    = "${meta.trimmer}-${meta.id}"
    maxmem = task.memory.toGiga()
    if ( params.spades_fix_cpus == -1 || task.cpus == params.spades_fix_cpus )
        """
        metaspades.py \
            ${params.spades_options} \
            --threads "${task.cpus}" \
            --only-assembler \
            --memory $maxmem \
            --pe1-1 ${reads[0]} \
            --pe1-2 ${reads[1]} \
            -o spades
        mv spades/assembly_graph_with_scaffolds.gfa ${prefix}_graph.gfa
        mv spades/scaffolds.fasta ${prefix}_scaffolds.fasta
        mv spades/contigs.fasta ${prefix}_contigs.fasta
        mv spades/spades.log ${prefix}.log
        gzip -c "${prefix}_contigs.fasta" > "${prefix}_contigs.fasta.gz"
        gzip "${prefix}_graph.gfa"
        gzip -c "${prefix}_scaffolds.fasta" > "${prefix}_scaffolds.fasta.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //g')
            metaspades: \$(metaspades.py --version | sed "s/SPAdes genome assembler v//; s/ \\[.*//")
        END_VERSIONS
        """
    else
        error "ERROR: '--spades_fix_cpus' was specified, but not succesfully applied. Likely this is caused by changed process properties in a custom config file."
}