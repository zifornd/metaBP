/*
 * Bowtie2 for read removal
 */
process BOWTIE2_REMOVAL_ALIGN {
    tag "${meta.trimmer}-${meta.id}-${task.ext.prefix}"
    
    conda (params.enable_conda ? "bioconda::bowtie2=2.4.2 bioconda::samtools=1.11 conda-forge::pigz=2.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' :
        'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.unmapped*.fastq.gz") , emit: reads
    path  "*.mapped*.read_ids.txt", optional:true , emit: read_ids
    tuple val(meta), path("*.bowtie2.log")        , emit: log
    path  'versions.yml'                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ? "${meta.trimmer}-${meta.id}-${task.ext.prefix}" : "${meta.trimmer}-${meta.id}"
    def sensitivity = params.host_removal_verysensitive ? "--very-sensitive" : "--sensitive"
    def save_ids = params.host_removal_save_ids ? "Y" : "N"
    if (!meta.single_end){
        """
        bowtie2 -p ${task.cpus} \
                -x ${index[0].getSimpleName()} \
                -1 "${reads[0]}" -2 "${reads[1]}" \
                $args \
                --un-conc-gz ${prefix}.unmapped_%.fastq.gz \
                --al-conc-gz ${prefix}.mapped_%.fastq.gz \
                1> /dev/null \
                2> ${prefix}.bowtie2.log
        if [ ${save_ids} = "Y" ] ; then
            gunzip -c ${prefix}.mapped_1.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.mapped_1.read_ids.txt
            gunzip -c ${prefix}.mapped_2.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.mapped_2.read_ids.txt
        fi
        rm -f ${prefix}.mapped_*.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        bowtie2 -p ${task.cpus} \
                -x ${index[0].getSimpleName()} \
                -U ${reads} \
                $args \
                --un-gz ${prefix}.unmapped.fastq.gz \
                --al-gz ${prefix}.mapped.fastq.gz \
                1> /dev/null \
                2> ${prefix}.bowtie2.log
        if [ ${save_ids} = "Y" ] ; then
            gunzip -c ${prefix}.mapped.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.mapped.read_ids.txt
        fi
        rm -f ${prefix}.mapped.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        END_VERSIONS
        """
    }
}