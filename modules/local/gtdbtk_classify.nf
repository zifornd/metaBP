// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process GTDBTK_CLASSIFY {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['assembler', 'id']) }

    conda (params.enable_conda ? "bioconda::gtdbtk=2.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gtdbtk:2.1.0--pyhdfd78af_5"
    } else {
        container "quay.io/biocontainers/gtdbtk:2.1.0--pyhdfd78af_5"
    }

    input:
    tuple val(meta), path("bins/*")
    tuple val(db_name), path("database/*")

    output:
    path "gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.*.summary.tsv"                 , emit: summary
    path "classify/gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.*.classify.tree.gz"   , emit: tree
    path "identify/gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.*.markers_summary.tsv", emit: markers
    path "align/gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.*.msa.fasta.gz"          , emit: msa
    path "align/gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.*.user_msa.fasta.gz"     , emit: user_msa
    path "align/gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.*.filtered.tsv"          , emit: filtered
    path "gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.log"                           , emit: log
    path "gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.warnings.log"                  , emit: warnings
    path "identify/gtdbtk.${meta.assembler}-${meta.trimmer}-${meta.id}.failed_genomes.tsv"   , emit: failed
    path '*.version.txt'                                                                     , emit: version

    script:
    def software = getSoftwareName(task.process)
    def pplacer_scratch = params.gtdbtk_pplacer_scratch ? "--scratch_dir pplacer_tmp" : ""
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
    """
    export GTDBTK_DATA_PATH="\${PWD}/database"
    if [ ${pplacer_scratch} != "" ] ; then
        mkdir pplacer_tmp
    fi

    gtdbtk classify_wf $options.args \
                    --genome_dir bins \
                    --prefix "gtdbtk.${prefix}" \
                    --out_dir "\${PWD}" \
                    --cpus ${task.cpus} \
                    --pplacer_cpus ${params.gtdbtk_pplacer_cpus} \
                    ${pplacer_scratch} \
                    --min_perc_aa ${params.gtdbtk_min_perc_aa} \
                    --min_af ${params.gtdbtk_min_af}

    gzip "classify/gtdbtk.${prefix}".*.classify.tree #> "gtdbtk.${prefix}".*.classify.tree.gz 
    #gzip "align/gtdbtk.${prefix}".*.msa.fasta > "gtdbtk.${prefix}".*.msa.fasta.gz
    mv gtdbtk.log "gtdbtk.${prefix}.log"
    mv gtdbtk.warnings.log "gtdbtk.${prefix}.warnings.log"
    gtdbtk --version | sed "s/gtdbtk: version //; s/ Copyright.*//" > ${software}.version.txt
    """
}
