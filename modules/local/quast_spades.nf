// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QUAST_SPADES {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['assembler', 'id']) }

    conda (params.enable_conda ? "bioconda::quast=5.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2"
    } else {
        container "quay.io/biocontainers/quast:5.0.2--py37pl526hb5aa323_2"
    }

    input:
    tuple val(meta), path(assembly)//, path(contigs_gz)
    tuple val(meta), path(contig)

    output:
    path "QUAST/*" , emit: qc
    path '*.version.txt'   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
       """
       metaquast.py --threads "${task.cpus}" --max-ref-number 0 -l "${prefix}_scaffolds","${prefix}_contigs" "${assembly}" "${contig}" -o "QUAST"
       metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//" > ${software}.version.txt
       """
}
