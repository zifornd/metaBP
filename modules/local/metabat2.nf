
process METABAT2 {
    tag "${meta.assembler}-${meta.trimmer}-${meta.id}"

    conda (params.enable_conda ? "bioconda::metabat2=2.15 conda-forge::python=3.6.7 conda-forge::biopython=1.74 conda-forge::pandas=1.1.5" : null)
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' :
        'quay.io/biocontainers/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' }"
   

    input:
    tuple val(meta), path(assembly), path(bam), path(bai)

    output:
    tuple val(meta), path("MetaBAT2/*.fa")                              , emit: bins
    path "${meta.assembler}-${meta.trimmer}-${meta.id}-depth.txt.gz"    , emit: depths
    path "MetaBAT2/discarded/*"                                         , emit: discarded
    path 'versions.yml'                                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.assembler}-${meta.trimmer}-${meta.id}"
    """
    OMP_NUM_THREADS=${task.cpus} jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
    metabat2 -t "${task.cpus}" -i "${assembly}" -a depth.txt -o "MetaBAT2/${prefix}" -m ${params.min_contig_size} --unbinned --seed ${params.metabat_rng_seed}

    gzip depth.txt
    mv depth.txt.gz "${prefix}-depth.txt.gz"

    # save unbinned contigs above thresholds into individual files, dump others in one file
    split_fasta.py "MetaBAT2/${prefix}.unbinned.fa" ${params.min_length_unbinned_contigs} ${params.max_unbinned_contigs} ${params.min_contig_size}

    # delete splitted file so that it doesnt end up in following processes
    rm "MetaBAT2/${prefix}.unbinned.fa"

    mkdir MetaBAT2/discarded
    gzip "MetaBAT2/${prefix}.lowDepth.fa" \
        "MetaBAT2/${prefix}.tooShort.fa" \
        "MetaBAT2/${prefix}.unbinned.pooled.fa" \
        "MetaBAT2/${prefix}.unbinned.remaining.fa"
    mv "MetaBAT2/${prefix}".*.fa.gz MetaBAT2/discarded/

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
END_VERSIONS
"""
}
