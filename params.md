# zifornd/metaBP pipeline parameters

Assembly, binning and annotation of metagenomes

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Input FastQ files or CSV samplesheet file containing information about the samples in the experiment. <details><summary>Help</summary><small>Use 
this to specify the location of your input FastQ files. For example:

```bash
--input 'path/to/data/sample_*_{1,2}.fastq.gz'
``` 

Alternatively, to assign different groups or to include long reads for hybrid assembly with metaSPAdes, you can specify a CSV samplesheet input file with 5 
columns and the following header: sample,group,short_reads_1,short_reads_2,long_reads. See 
(https://nf-co.re/mag/usage#input-specifications).</small></details>| `string` |  |  |  |
| `single_end` | Specifies that the input is single-end reads. <details><summary>Help</summary><small>By default, the pipeline expects paired-end data. If 
you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation 
marks, can then be used for `--input`. For example:

```bash
--single_end --input '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.</small></details>| `boolean` |  |  |  |
| `outdir` | Path to the output directory where the results will be saved. | `string` | . |  |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail 
with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on 
the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `igenomes_base` | Directory / URL base for iGenomes references. | `string` | s3://ngi-igenomes/igenomes |  | True |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the 
pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>| `boolean` |
|  | True |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be 
able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download 
the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | 
https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `hostnames` | Institutional configs hostname. | `string` |  |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU 
requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the 
memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 62.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time
requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option 
specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. 
See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary
email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `tracedir` | Directory to keep pipeline Nextflow logs and reports. | `string` | ${params.outdir}/pipeline_info |  | True |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `show_hidden_params` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are 
not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| 
`boolean` |  |  | True |
| `enable_conda` | Run this workflow with Conda. You can also use '-profile conda' instead of providing this parameter. | `boolean` |  |  | True |
| `singularity_pull_docker_container` | Instead of directly downloading Singularity images for use with Singularity, force the workflow to pull and convert 
Docker containers instead. <details><summary>Help</summary><small>This may be useful for example if you are unable to directly pull Singularity containers to
run the pipeline due to http/https proxy issues.</small></details>| `boolean` |  |  | True |

## Reproducibility options

Use these parameters to also enable reproducible results from the individual assembly and binning tools .

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `megahit_fix_cpu_1` | Fix number of CPUs for MEGAHIT to 1. Not increased with retries. <details><summary>Help</summary><small>MEGAHIT only generates 
reproducible results when run single-threaded. 

When using this parameter do not change the number of CPUs for the `megahit` process with a custom config file. This would result in an error.

Default: The number of CPUs is specified in the `base.config` file, and increased with each retry.</small></details>| `boolean` |  |  |  |
| `spades_fix_cpus` | Fix number of CPUs used by SPAdes. Not increased with retries. <details><summary>Help</summary><small>SPAdes is designed to be 
deterministic for a given number of threads. To generate reproducible results fix the number of CPUs using this parameter.

When using this parameter do not change the number of CPUs for the `spades` process with a custom config file. This would result in an error.

Default: -1 (the number of CPUs is specified in the `base.config` or in a custom config file, and increased with each retry).</small></details>| `integer` | 
-1 |  |  |
| `metabat_rng_seed` | RNG seed for MetaBAT2. <details><summary>Help</summary><small>MetaBAT2 is run by default with a fixed seed within this pipeline, thus 
producing reproducible results. You can set it also to any other positive integer to ensure reproducibility. Set the parameter to 0 to use a random 
seed.</small></details>| `integer` | 1 |  |  |

## Quality control for short reads options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `cutadapt` | Run Cutadapt for quality trimming. | `boolean` |  |  |  |
| `cutadapt_quality_cutoff` |  | `integer` | 20 |  |  |
| `cutadapt_min_read_length` |  | `integer` | 50 |  |  |
| `cutadapt_max_error_rate` |  | `number` | 0.25 |  |  |
| `trimmomatic` | Run Trimmomatic for quality trimming. | `boolean` | True |  |  |
| `trimmomatic_illuminaclip` |  | `string` | NexteraPE-PE.fa:2:30:10:8:TRUE |  |  |
| `trimmomatic_sliding_window` |  | `string` | 4:20 |  |  |
| `trimmomatic_min_length` |  | `integer` | 50 |  |  |
| `save_trimmed_fail` | Save the trimmed FastQ files in the results directory. <details><summary>Help</summary><small>By default, trimmed FastQ files will 
not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when 
complete.</small></details>| `boolean` |  |  |  |
| `host_genome` | Name of iGenomes reference for host contamination removal. <details><summary>Help</summary><small>This parameter is mutually exclusive with
`--host_genome`. Host read removal is done with Bowtie2. 
Both the iGenomes FASTA file as well as corresponding, already pre-built Bowtie 2 index files will be used.</small></details>| `string` |  |  |  |
| `host_fasta` | Fasta reference file for host contamination removal. <details><summary>Help</summary><small>This parameter is mutually exclusive with 
`--host_fasta`. The reference can be masked. Host read removal is done with Bowtie2.</small></details>| `string` |  |  |  |
| `host_removal_verysensitive` | Use the `--very-sensitive` instead of the`--sensitive`setting for Bowtie 2 to map reads against the host genome. | `boolean`
|  |  |  |
| `host_removal_save_ids` | Save the read IDs of removed host reads. | `boolean` |  |  |  |
| `keep_phix` | Keep reads similar to the Illumina internal standard PhiX genome. | `boolean` |  |  |  |
| `phix_reference` | Genome reference used to remove Illumina PhiX contaminant reads. | `string` | 
${baseDir}/assets/data/GCA_002596845.1_ASM259684v1_genomic.fna.gz |  | True |

## Taxonomic profiling options

Taxonomic classification is disabled by default. You have to specify one of the options below to activate it.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `kraken2` | Run Kraken2 for taxonomic classification of reads. | `boolean` | True |  |  |
| `kraken2_db` | Database for taxonomic binning with kraken2. <details><summary>Help</summary><small>The database file must be a compressed tar archive that 
contains at least the three files `hash.k2d`, `opts.k2d` and `taxo.k2d`. E.g. 
ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz.</small></details>| `string` |  |  |  |
| `skip_krona` | Skip creating a krona plot for taxonomic binning. | `boolean` |  |  |  |
| `gtdbtk` | Run GTDB-Tk for taxonomic classification of bins. | `boolean` |  |  |  |
| `gtdb` | GTDB database for taxonomic classification of bins with GTDB-tk. <details><summary>Help</summary><small>For information which GTDB reference 
databases are compatible with the used GTDB-tk version see 
https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data.</small></details>| `string` | 
https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz |  |  |
| `gtdbtk_min_completeness` | Min. bin completeness (in %) required to apply GTDB-tk classification. <details><summary>Help</summary><small>Completeness 
assessed with BUSCO analysis (100% - %Missing). Must be greater than 0 (min. 0.01) to avoid GTDB-tk errors. If too low, GTDB-tk classification results can be
impaired due to not enough marker genes!</small></details>| `number` | 50 |  |  |
| `gtdbtk_max_contamination` | Max. bin contamination (in %) allowed to apply GTDB-tk classification. <details><summary>Help</summary><small>Contamination 
approximated based on BUSCO analysis (%Complete and duplicated). If too high, GTDB-tk classification results can be impaired due to 
contamination!</small></details>| `number` | 10 |  |  |
| `gtdbtk_min_perc_aa` | Min. fraction of AA (in %) in the MSA for bins to be kept. | `number` | 10 |  |  |
| `gtdbtk_min_af` | Min. alignment fraction to consider closest genome. | `number` | 0.65 |  |  |
| `gtdbtk_pplacer_cpus` | Number of CPUs used for the by GTDB-Tk run tool pplacer. <details><summary>Help</summary><small>A low number of CPUs helps to 
reduce the memory required/reported by GTDB-Tk. See also the [GTDB-Tk 
documentation](https://ecogenomics.github.io/GTDBTk/faq.html#gtdb-tk-reaches-the-memory-limit-pplacer-crashes).</small></details>| `number` | 1 |  |  |
| `gtdbtk_pplacer_scratch` | Reduce GTDB-Tk memory consumption by running pplacer in a setting writing to disk. <details><summary>Help</summary><small>Will 
be slower. Set to `false` to turn this off.</small></details>| `boolean` | True |  |  |
| `sourmash` | Run Sourmash or taxonomic classification of bins. | `boolean` | True |  |  |
| `sourmash_db` | Sourmash LCA database. <details><summary>Help</summary><small>The database file must be an LCA(.lca.json.gz) database.</small></details>| 
`string` |  |  |  |
| `sourmash_sketch_scaled` | Create a scaled MinHash with k-mers sampled deterministically at 1 per <scaled> value. | `number` | 1000 |  |  |
| `sourmash_sketch_k` | k-mer size to create a sketch at this size. <details><summary>Help</summary><small>Typically ksize is between 4 and 
100.</small></details>| `number` | 31 |  |  |

## Assembly options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `coassemble_group` | Co-assemble samples within one group, instead of assembling each sample separately. | `boolean` |  |  |  |
| `spades_options` | Additional custom options for SPAdes. <details><summary>Help</summary><small>An example is adjusting k-mers ("-k 21,33,55,77") or adding
(https://github.com/ablab/spades#advanced-options). But not -t, -m, -o or --out-prefix, because these are already in use.</small></details>| `string` |  |  |
|
| `megahit_options` | Additional custom options for MEGAHIT. <details><summary>Help</summary><small>An example is adjusting presets (e.g. "--presets 
meta-large"), k-mers (e.g. "-k 21,33,55,77") or adding other (https://github.com/voutcn/megahit#advanced-usage). For example, increase the minimum k-mer in 
the event of an error message such as "Too many vertices in the unitig graph, you may increase the kmer size to remove tons of erroneous kmers." in the 
MEGAHIT log file. But not --threads, --memory, -o or input read files, because these are already in use.</small></details>| `string` |  |  |  |
| `skip_spades` | Skip Illumina-only SPAdes assembly. | `boolean` |  |  |  |
| `skip_megahit` | Skip MEGAHIT assembly. | `boolean` |  |  |  |
| `skip_quast` | Skip metaQUAST. | `boolean` |  |  |  |

## Gene prediction options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_prodigal` | Skip Prodigal gene prediction | `boolean` |  |  |  |

## Binning options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `binning_map_mode` | Defines mapping strategy to compute co-abundances for binning, i.e. which samples will be mapped against the assembly. 
<details><summary>Help</summary><small>Available: `all`, `group` or `own`. Note that `own` cannot be specified in combination with `--coassemble_group`.

Note that specifying `all` without additionally specifying `--coassemble_group` results in `n^2` mapping processes for each assembly method, where `n` is the
number of samples.</small></details>| `string` | group |  |  |
| `skip_binning` | Skip metagenome binning. | `boolean` |  |  |  |
| `min_contig_size` | Minimum contig size to be considered for binning and for bin quality check. <details><summary>Help</summary><small>For forwarding into 
downstream analysis, i.e. QUAST and BUSCO, and reporting.</small></details>| `integer` | 1500 |  |  |
| `min_length_unbinned_contigs` | Minimal length of contigs that are not part of any bin but treated as individual genome. 
<details><summary>Help</summary><small>Contigs that do not fulfill the thresholds of `--min_length_unbinned_contigs` and `--max_unbinned_contigs` are pooled 
for downstream analysis and reporting, except contigs that also do not fullfill `--min_contig_size` are not considered further.</small></details>| `integer` 
| 1000000 |  |  |
| `max_unbinned_contigs` | Maximal number of contigs that are not part of any bin but treated as individual genome. 
<details><summary>Help</summary><small>Contigs that do not fulfill the thresholds of `--min_length_unbinned_contigs` and `--max_unbinned_contigs` are pooled 
for downstream analysis and reporting, except contigs that also do not fullfill `--min_contig_size` are not considered further.</small></details>| `integer` 
| 100 |  |  |
| `skip_prokka` | Skip Prokka genome annotation. | `boolean` |  |  |  |

## Bin quality check options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_busco` | Disable bin QC with BUSCO. | `boolean` |  |  |  |
| `busco_reference` | Download path for BUSCO lineage dataset, instead of using automated lineage selection. <details><summary>Help</summary><small>E.g. 
https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2020-03-06.tar.gz. Available databases are listed here: 
https://busco-data.ezlab.org/v5/data/lineages/.</small></details>| `string` |  |  |  |
| `busco_download_path` | Path to local folder containing already downloaded and unpacked lineage datasets. <details><summary>Help</summary><small>If 
provided, BUSCO analysis will be run in offline mode. Data can be downloaded from https://busco-data.ezlab.org/v5/data/ (files still need to be unpacked 
manually). Run in combination with automated lineage selection.</small></details>| `string` |  |  |  |
| `busco_auto_lineage_prok` | Run BUSCO with automated lineage selection but ignoring eukaryotes (saves runtime). | `boolean` |  |  |  |
| `save_busco_reference` | Save the used BUSCO lineage datasets provided via --busco_reference or downloaded when not using --busco_reference or 
--busco_download_path. <details><summary>Help</summary><small>Useful to allow reproducibility, as BUSCO datasets are frequently updated and old versions do 
not always remain accessible.</small></details>| `boolean` |  |  |  |

