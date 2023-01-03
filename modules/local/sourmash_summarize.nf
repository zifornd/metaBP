/*
--------------------------Sourmash added for Taxonomic Classification by Zifo----------------------
*/


// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process SOURMASH_SUMMARIZE {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"

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
    path(db)
    tuple val(meta), path(signatures)

    output:
    tuple val(meta), path('*.csv')	    , emit: report
    path "*.version.txt"                    , emit: version
    
    script:
    def software = getSoftwareName(task.process)
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
    
	"""
    sourmash lca summarize \\
        $options.args \\
        --db ${db} \\
	--query \\
	${signatures} \\
        -o '${prefix}.csv' 
    	echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' > ${software}.version.txt
        """
}

