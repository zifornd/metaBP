/*
----------------------Trimmomatic added for Quality Control by Zifo----------------------
*/
process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'
    
    /* conda (params.enable_conda ? 'bioconda::trimmomatic=0.36' : null)
    if ( workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container){
        container 'https://depot.galaxyproject.org/singularity/trimmomatic:0.32--hdfd78af_4'
    } else {
       container 'quay.io/biocontainers/trimmomatic:0.38--0'
    } */
    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path "versions.yml"                     , emit: versions
    //tuple val(meta), path('*.fail.fastq.gz'), optional:true, emit: reads_fail

    when:
    task.ext.when == null || task.ext.when

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}":"${meta.id}"
    def type = meta.single_end ? "SE" : "PE"
    def trimmed  = meta.single_end ? "${prefix}.trim.fastq.gz" 
        : "${prefix}_1.trim.fastq.gz ${prefix}_1.single.trim.1.fastq.gz ${prefix}_2.trim.fastq.gz ${prefix}_2.single.trim.2.fastq.gz"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    trimmomatic $type \\
        -threads $task.cpus \\
        -phred33 \\
        ${reads} \\
        ${trimmed} \\
        $args \\
        2> ${prefix}.trimmomatic.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}

