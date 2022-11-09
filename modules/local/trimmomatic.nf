/*
----------------------Trimmomatic added for Quality Control by Zifo----------------------
*/


// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::trimmomatic=0.36' : null)
    if ( workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container){
        container 'https://depot.galaxyproject.org/singularity/trimmomatic:0.32--hdfd78af_4'
    } else {
       container 'quay.io/biocontainers/trimmomatic:0.38--0'
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path "*.version.txt"                    , emit: version
    tuple val(meta), path('*.fail.fastq.gz'), optional:true, emit: reads_fail


    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}":"${meta.id}"
    if (meta.single_end) {
        def fail_fastq = params.save_trimmed_fail ? "--failed_out ${prefix}.fail.fastq.gz" : ''
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz

        trimmomatic SE \\
        -threads $task.cpus \\
	-phred33 \\
        ${prefix}.fastq.gz \\
	${prefix}.trim.fastq.gz \\
        $fail_fastq \\
	$options.args \\
        2> ${prefix}.trimmomatic.log
	trimmomatic -version > ${software}.version.txt
        #trimmomatic -version 2>&1 | sed -e "s/trimmomatic //g" > ${software}.version.txt
        """

      } else {
        def fail_fastq = params.save_trimmed_fail ? "--unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trimmomatic PE \\
        -threads $task.cpus \\
        -phred33 \\
	${prefix}_1.fastq.gz \\
        ${prefix}_2.fastq.gz \\
        ${prefix}_1.trim.fastq.gz \\
	${prefix}_1.single.trim.1.fastq.gz \\
        ${prefix}_2.trim.fastq.gz \\
	${prefix}_2.single.trim.2.fastq.gz \\
        $options.args \\
	$fail_fastq \\
        2> ${prefix}.log
	trimmomatic -version > ${software}.version.txt
        #echo \$(trimmomatic -version 2>&1) | sed -e "s/trimmomatic //g" > ${software}.version.txt
        """
       }
    }

