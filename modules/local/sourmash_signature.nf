/*
--------------------------Sourmash added for Taxonomic Classification by Zifo----------------------
*/


// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process SOURMASH_SIGNATURE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::sourmash=4.4.0' : null)
    if ( workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container){
        container 'https://depot.galaxyproject.org/singularity/sourmash:4.4.0--hdfd78af_0'
    } else {
       container 'quay.io/biocontainers/sourmash:4.4.0--hdfd78af_0'
    }
	

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path('*.sig')	    , emit: signatures
    //tuple val(meta), path('*.log')          , emit: log
    path "*.version.txt"                    , emit: version
    
    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.assembler}-${meta.id}-${options.suffix}":"${meta.id}"
    
	"""
        sourmash sketch dna ${options.args} ${bins} --output "${prefix}.sig"      
    	sourmash --version | sed "s/SOURMASH //" > ${software}.version.txt
   	"""
}

