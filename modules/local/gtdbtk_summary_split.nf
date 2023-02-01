/*
--------------------------Module added to split GTDB-Tk summary file by Zifo----------------------
*/


// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process GTDBTK_SUMMARY_SPLIT {
    //tag "${meta.assembler}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['assembler', 'id']) }

    input:
    file(gtdbtk_sum)

    output:
    path("spades_cutadapt_bin_summary.tsv")    , emit: cut_spa_gtdbsum
    path("megahit_cutadapt_bin_summary.tsv")   , emit: cut_mg_gtdbsum
    path("spades_trimmomatic_bin_summary.tsv") , emit: trim_spa_gtdbsum
    path("megahit_trimmomatic_bin_summary.tsv"), emit: trim_mg_gtdbsum

    script:
    //def software = getSoftwareName(task.process)
    """
    grep -i 'user_genome\\|cutadapt' ${gtdbtk_sum} > cutadapt_bin_summary.tsv 
    grep -i 'user_genome\\|spades' cutadapt_bin_summary.tsv > spades_cutadapt_bin_summary.tsv
    grep -i 'user_genome\\|megahit' cutadapt_bin_summary.tsv > megahit_cutadapt_bin_summary.tsv

    grep -i 'user_genome\\|trimmomatic' ${gtdbtk_sum} > trimmomatic_bin_summary.tsv
    grep -i 'user_genome\\|spades' trimmomatic_bin_summary.tsv > spades_trimmomatic_bin_summary.tsv
    grep -i 'user_genome\\|megahit' trimmomatic_bin_summary.tsv > megahit_trimmomatic_bin_summary.tsv
    """
 }