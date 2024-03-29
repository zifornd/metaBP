# MetaBP: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the Intermediate_Outputs and Results directory after the pipeline has finished. The tools giving the final output of the pipeline are stored in Results directory whereas rest other tools giving intermediate outputs are stored in Intermediate_Outputs directory. All paths are relative to the top-level Intermediate_Outputs and Results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Quality control](#quality-control) of input reads - trimming and contaminant removal
* [Taxonomic classification of trimmed reads](#taxonomic-classification-of-trimmed-reads)
* [Assembly](#assembly) of trimmed reads
* [Protein-coding gene prediction](#gene-prediction) of assemblies
* [Binning](#binning) of assembled contigs
* [Taxonomic classification of binned genomes](#taxonomic-classification-of-binned-genomes)
* [Genome annotation of binned genomes](#genome-annotation-of-binned-genomes)
* [Additional summary for binned genomes](#additional-summary-for-binned-genomes)
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

Note that when specifying the parameter `--coassemble_group`, for the corresponding output filenames/directories of the assembly or downsteam processes the group ID, or more precisely the term `group-[group_id]`, will be used instead of the sample ID.

## Quality control

These steps trim away the adapter sequences present in input reads, trims away bad quality bases and discard reads that are too short.
It also removes host contaminants and sequencing controls, such as PhiX or the Lambda phage.
FastQC is run for visualising the general quality metrics of the sequencing runs before and after trimming.

<!-- TODO: Add example MultiQC plots generated for this pipeline -->

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/FastQC/`
    * `[sample]_[1/2]_fastqc.html`: FastQC report, containing quality metrics for your untrimmed raw fastq files
    * `[sample]_[trimmer]_[1/2]_fastqc.html`: FastQC report, containing quality metrics for trimmed and, if specified, filtered read.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### fastp

[fastp](https://github.com/OpenGene/fastp) is a all-in-one fastq preprocessor for read/adapter trimming and quality control. It is used in this pipeline for trimming adapter sequences and discard low-quality reads. Its output is in the Intermediate_Outputs folder and part of the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/fastp/[sample/group]/`
    * `[sample].fastp.html`: Interactive report
    * `[sample].fastp.json`: Report in json format

</details>

### Cutadapt

[Cutadapt](https://github.com/marcelm/cutadapt/) finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from high-throughput sequencing reads. It is used in this pipeline for trimming adapter sequences and discard low-quality reads. Its output is in the Intermediate_Outputs folder and part of the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/Cutadapt//[sample/group]/`
    * `[sample].cutadapt.log`: 


</details>

### Trimmomatic

[Trimmomatic](https://github.com/usadellab/Trimmomatic) performs a variety of useful trimming tasks for paired-end and single ended data. It is used in this pipeline for trimming adapter sequences and discard low-quality reads. Its output is in the Intermediate_Outputs folder and part of the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/Trimmomatic/[sample/group]/`
    * `[sample].trimmomatic.log`: Contains a brief log file indicating the number of input reads, surviving reads and dropped reads.

</details>

### Remove PhiX sequences from short reads

The pipeline uses bowtie2 to map the reads against PhiX and removes mapped reads.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/phiX Removal/`
    * `[trimmer]-[sample]-phix_removed.bowtie2.log`: Contains a brief log file indicating how many reads have been retained.

</details>

### Host read removal

The pipeline uses bowtie2 to map short reads against the host reference genome specified with `--host_genome` or `--host_fasta` and removes mapped reads. The information about discarded and retained reads is also included in the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/Host Removal/` 
    * `[trimmer]-[sample]-host_removed.bowtie2.log`: Contains the bowtie2 log file indicating how many reads have been mapped as well as a file listing the read ids of discarded reads.

</details>



## Taxonomic classification of trimmed reads

### Kraken

[Kraken2](https://github.com/DerrickWood/kraken2) classifies reads using a k-mer based approach as well as assigns taxonomy using a Lowest Common Ancestor (LCA) algorithm.

<details markdown="1">
<summary>Output files</summary>

* `Results/Kraken2/[trimmer]/[sample]/`
    * `[trimmer]-[sample]_kraken2.report.txt`: Classification in the Kraken report format. See the [kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) for more details
     
    
</details>

### Krona

[Krona](https://github.com/marbl/Krona) allows hierarchical data to be explored with zooming, multi-layered pie charts. Krona charts can be created using an Excel template or KronaTools, which includes support for several bioinformatics tools and raw data formats. The interactive charts are self-contained and can be viewed with any modern web browser.

<details markdown="1">
<summary>Output files</summary>

* `Results/Krona/[trimmer]/kraken2/[sample]/`
    * `taxonomy.krona.html`: Interactive pie chart of the cutadapt output produced by [KronaTools](https://github.com/marbl/Krona/wiki)
     
</details>




## Assembly

Trimmed (short) reads are assembled with both megahit and SPAdes. Hybrid assembly is only supported by SPAdes.

### MEGAHIT

[MEGAHIT](https://github.com/voutcn/megahit) is a single node assembler for large and complex metagenomics short reads.

<details markdown="1">
<summary>Output files</summary>


* `Intermediate_Outputs/MEGAHIT/[trimmer]/`
    * `MEGAHIT/[trimmer]-[sample/group].contigs.fa.gz`: Compressed metagenome assembly in fasta format
    * `MEGAHIT/[trimmer]-[sample/group].log`: Log file
    
    * `QC/[sample/group]/`: Directory containing QUAST files and Bowtie2 mapping logs
        * `MEGAHIT-[sample/group]-[trimmer].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
        * `MEGAHIT-[sample/group]-[sampleToMap]-[trimmer].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").

</details>

### SPAdes

[SPAdes](http://cab.spbu.ru/software/spades/) was originally a single genome assembler that later added support for assembling metagenomes.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/SPAdes/[trimmer]/`
    * `[trimmer]-[sample/group]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
    * `[trimmer]-[sample/group]_graph.gfa.gz`: Compressed assembly graph in gfa format
    * `[trimmer]-[sample/group]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
    * `[trimmer]-[sample/group]pt.log`: Log file
    * `QC/[sample/group]/`: Directory containing QUAST files and Bowtie2 mapping logs
        * `SPAdes-[sample/group]-[trimmer].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the trimmomatic output sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
        * `SPAdes-[sample/group]-[sampleToMap]-[trimmer].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").

</details>


### Metagenome QC with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates metagenome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/QUAST/[assembler]-[trimmer]-[sample/group].bin/`
    * `report.*`: QUAST report in various formats, such as html, txt, tsv or tex
    * `quast.log`: QUAST log file
    * `predicted_genes/[assembler]-[trimmer]-[sample/group].rna.gff`: Contig positions for rRNA genes in gff version 3 format

</details>

## Gene prediction

Protein-coding genes are predicted for each assembly.

### Prodigal

[Prodigal](https://github.com/hyattpd/Prodigal) is a protein-coding gene prediction software tool for bacterial and archaeal genomes. 

<details markdown="1">
<summary>Output files</summary>

* `Results/Prodigal/[assembler]/[sample/group]/`
    * `[sample/group].gff`: Gene Coordinates in GFF format
    * `[sample/group].faa`: The protein translation file consists of all the proteins from all the sequences in multiple FASTA format.
    * `[sample/group].fna`: Nucleotide sequences of the predicted proteins using the DNA alphabet, not mRNA (so you will see 'T' in the output and not 'U').
    * `[sample/group]_all.txt`: Information about start positions of genes.

</details>




## Binning

### Contig sequencing depth

Sequencing depth per contig and sample is generated by `jgi_summarize_bam_contig_depths --outputDepth`. The values correspond to `(sum of exactely aligned bases) / ((contig length)-2*75)`. For example, for two reads aligned exactly with `10` and `9` bases on a 1000 bp long contig the depth is calculated by `(10+9)/(1000-2*75)` (1000bp length of contig minus 75bp from each end, which is excluded).

<details markdown="1">
<summary>Output files</summary>

* `Binning/[trimmer]/`
    * `[assembler]-[trimmer]-[sample/group]-depth.txt.gz`: Sequencing depth for each contig and sample or group, only for short reads.

</details>

### MetaBAT2

[MetaBAT2](https://bitbucket.org/berkeleylab/metabat) recovers genome bins (that is, contigs/scaffolds that all belongs to a same organism) from metagenome assemblies.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/Binning/[trimmer]/MetaBAT2/`
    * `[assembler]-[trimmer]-[sample/group].*.fa`: Genome bins retrieved from input assembly
    * `[assembler]-[trimmer]-[sample/group].unbinned.*.fa`: Contigs that were not binned with other contigs but considered interesting. By default, these are at least 1 Mbp (`--min_length_unbinned_contigs`) in length and at most the 100 longest contigs (`--max_unbinned_contigs`) are reported

</details>

All the files and contigs in this folder will be assessed by QUAST and BUSCO.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/Binning/[trimmer]/MetaBAT2/discarded/`
    * `*.lowDepth.fa.gz`: Low depth contigs that are filtered by MetaBat2
    * `*.tooShort.fa.gz`: Too short contigs that are filtered by MetaBat2
    * `*.unbinned.pooled.fa.gz`: Pooled unbinned contigs equal or above `--min_contig_size`, by default 1500 bp.
    * `*.unbinned.remaining.fa.gz`: Remaining unbinned contigs below `--min_contig_size`, by default 1500 bp, but not in any other file.

</details>

All the files in this folder contain small and/or unbinned contigs that are not further processed.

Files in these two folders contain all contigs of an assembly.

### Bin sequencing depth

For each genome bin the median sequencing depth is computed based on the corresponding contig depths given in `Binning/[assembler]-[trimmer]-[sample/group]-depth.txt.gz`.

<details markdown="1">
<summary>Output files</summary>

* `Binning/[trimmer]/`
    * `bin_depths_summary.tsv`: Summary of bin sequencing depths for all samples. Depths are available for samples mapped against the corresponding assembly, i.e. according to the mapping strategy specified with `--binning_map_mode`. Only for short reads.
    * `[assembler]-[sample/group]-binDepths.heatmap.png`: Clustered heatmap showing bin abundances of the assembly across samples. Bin depths are transformed to centered log-ratios and bins as well as samples are clustered by Euclidean distance. Again, sample depths are available according to the mapping strategy specified with `--binning_map_mode`.

</details>

### QC for metagenome assembled genomes with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates genome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/QUAST/[assembler]-[trimmer]-[sample/group].bin/`
    * `report.*`: QUAST report in various formats, such as html, txt, tsv or tex
    * `quast.log`: QUAST log file
    * `predicted_genes/[assembler]-[sample/group].rna.gff`: Contig positions for rRNA genes in gff version 3 format
* `Intermediate_Outputs/QUAST/[assembler]-[trimmer]-[sample/group]-quast_summary.tsv`: QUAST output for all bins summarized

</details>

### QC for metagenome assembled genomes with BUSCO

[BUSCO](https://busco.ezlab.org/) is a tool used to assess the completeness of a genome assembly. It is run on all the genome bins and high quality contigs obtained by MetaBAT2. By default, BUSCO is run in automated lineage selection mode in which it first tries to select the domain and then a more specific lineage based on phylogenetic placement. If available, result files for both the selected domain lineage and the selected more specific lineage are placed in the output directory. If a lineage dataset is specified already with `--busco_reference`, only results for this specific lineage will be generated.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/BUSCO/[assembler]/`
    * `[assembler]-[trimmer]-[bin]_busco.log`: Log file containing the standard output of BUSCO.
    * `[assembler]-[trimmer]-[bin]_busco.err`: File containing potential error messages returned from BUSCO.
    * `short_summary.domain.[lineage].[assembler]-[bin].txt`: BUSCO summary of the results for the selected domain when run in automated lineage selection mode. Not available for bins for which a viral lineage was selected.
    * `short_summary.specific_lineage.[lineage].[assembler]-[bin].txt`: BUSCO summary of the results in case a more specific lineage than the domain could be selected or for the lineage provided via `--busco_reference`.
    * `[assembler]-[trimmer]-[bin]_buscos.[lineage].fna.gz`: Nucleotide sequence of all identified BUSCOs for used lineages (domain or specific).
    * `[assembler]-[trimmer]-[bin]_buscos.[lineage].faa.gz`: Aminoacid sequence of all identified BUSCOs for used lineages (domain or specific).
    * `[assembler]-[trimmer]-[bin]_prodigal.gff`: Genes predicted with Prodigal.

</details>

If the parameter `--save_busco_reference` is set, additionally the used BUSCO lineage datasets are stored in the output directy.

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/BUSCO/`
    * `busco_downloads/`: All files and lineage datasets downloaded by BUSCO when run in automated lineage selection mode. (Can currently not be used to reproduce analysis.
    * `reference/*.tar.gz`: BUSCO reference lineage dataset that was provided via `--busco_reference`.

</details>

Besides the reference files or output files created by BUSCO, the following summary files will be generated:

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/BUSCO/`
    * `[assembler]-[trimmer]-busco_summary.tsv`: A summary table of the BUSCO results, with % of marker genes found. If run in automated lineage selection mode, both the results for the selected domain and for the selected more specific lineage will be given, if available.

</details>

## Taxonomic classification of binned genomes


### Sourmash

[Sourmash](https://github.com/sourmash-bio/sourmash) is a k-mer based taxonomic exploration and classification routines for genome and metagenome analysis.

<details markdown="1">
<summary>Output files</summary>

* `Results/Sourmash/`
    * `[assembler]-[trimmer]-[sample/group].csv`: A csv file containing the classification of the bins.

</details>

### GTDB-Tk

[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) is a toolkit for assigning taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy [GTDB](https://gtdb.ecogenomic.org/). MetaBP uses GTDB-Tk to classify binned genomes which satisfy certain quality criteria (i.e. completeness and contamination assessed with the BUSCO analysis).

<details markdown="1">
<summary>Output files</summary>

* `Results/GTDB-Tk/[assembler]/`
    * `gtdbtk.[assembler]-[trimmer]-[sample/group].{bac120/ar122}.summary.tsv`: Classifications for bacterial and archaeal genomes (see the [GTDB-Tk documentation for details](https://ecogenomics.github.io/GTDBTk/files/summary.tsv.html).
    * `gtdbtk.[assembler]-[trimmer]-[sample/group].{bac120/ar122}.classify.tree.gz`: Reference tree in Newick format containing query genomes placed with pplacer.
    * `gtdbtk.[assembler]-[trimmer]-[sample/group].{bac120/ar122}.markers_summary.tsv`: A summary of unique, duplicated, and missing markers within the 120 bacterial marker set, or the 122 archaeal marker set for each submitted genome.
    * `gtdbtk.[assembler]-[trimmer]-[sample/group].{bac120/ar122}.msa.fasta.gz`: FASTA file containing MSA of submitted and reference genomes.
    * `gtdbtk.[assembler]-[trimmer]-[sample/group].{bac120/ar122}.filtered.tsv`: A list of genomes with an insufficient number of amino acids in MSA.
    * `gtdbtk.[assembler]-[trimmer]-[sample/group].*.log`: Log files.
    * `gtdbtk.[assembler]-[trimmer]-[sample/group].failed_genomes.tsv`: A list of genomes for which the GTDB-Tk analysis failed, e.g. because Prodigal could not detect any genes.
* `Results/GTDB-Tk/gtdbtk_summary-[assembler]-[trimmer].tsv`: A summary table of the GTDB-Tk classification results for all bins, also containing bins which were discarded based on the BUSCO QC, which were filtered out by GTDB-Tk ((listed in `*.filtered.tsv`) or for which the analysis failed (listed in `*.failed_genomes.tsv`).

</details>

## Genome annotation of binned genomes

### Prokka

Whole genome annotation is the process of identifying features of interest in a set of genomic DNA sequences, and labelling them with useful information. [Prokka](https://github.com/tseemann/prokka) is a software tool to annotate bacterial, archaeal and viral genomes quickly and produce standards-compliant output files.

<details markdown="1">
<summary>Output files</summary>

* `Results/Prokka/[assembler]/[assembler]-[trimmer]-[sample/group].[bin]/`
    * `[assembler]-[trimmer]-[bin].gff`: annotation in GFF3 format, containing both sequences and annotations
    * `[assembler]-[trimmer]-[bin].gbk`: annotation in GenBank format, containing both sequences and annotations
    * `[assembler]-[trimmer]-[bin].fna`: nucleotide FASTA file of the input contig sequences
    * `[assembler]-[trimmer]-[bin].faa`: protein FASTA file of the translated CDS sequences
    * `[assembler]-[trimmer]-[bin].ffn`: nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
    * `[assembler]-[trimmer]-[bin].sqn`: an ASN1 format "Sequin" file for submission to Genbank
    * `[assembler]-[trimmer]-[bin].fsa`: nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file
    * `[assembler]-[trimmer]-[bin].tbl`: feature Table file, used by "tbl2asn" to create the .sqn file
    * `[assembler]-[trimmer]-[bin].err`: unacceptable annotations - the NCBI discrepancy report.
    * `[assembler]-[trimmer]-[bin].log`: contains all the output that Prokka produced during its run
    * `[assembler]-[trimmer]-[bin].txt`: statistics relating to the annotated features found
    * `[assembler]-[trimmer]-[bin].tsv`: tab-separated file of all features (locus_tag, ftype, len_bp, gene, EC_number, COG, product)

</details>

## Additional summary for binned genomes

<details markdown="1">
<summary>Output files</summary>

* `Intermediate_Outputs/Binning/[assembler]-[trimmer]-[sample/group]/bin_summary_[trimmer]_[assembler].tsv`: Summary of bin sequencing depths together with BUSCO, QUAST and GTDB-Tk results, if at least one of the later was generated.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.tsv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
