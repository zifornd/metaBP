/*
--------------------------Cutadapt added for Quality Control by Zifo----------------------
*/
process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::cutadapt=3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
        'quay.io/biocontainers/cutadapt:3.4--py39h38f01e4_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path "versions.yml"                    , emit: versions
    //tuple val(meta), path('*.fail.fastq.gz'), optional:true, emit: reads_fail

    when:
    task.ext.when == null || task.ext.when

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}":"${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    cutadapt \\
        --cores $task.cpus \\
        $args \\
        $trimmed \\
        ${reads} \\
        > ${prefix}.cutadapt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}