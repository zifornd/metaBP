/*
============================ MetaBP - Shotgun Metagenomics Sequencing Tools Benchmarking Pipeline ==============================
  Written by Zifo RnD Solutions
  nf-core/mag pipeline 2.1.1 was modified 



/*
========================================================================================
                                   VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check already if long reads are provided
//include { hasExtension } from '../modules/local/functions'

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}
def hybrid = false
if(hasExtension(params.input, "csv")){
    Channel
        .from(file(params.input))
        .splitCsv(header: true)
        .map { row ->
                if (row.size() == 5) {
                    if (row.long_reads) hybrid = true
                } else {
                    log.error "Input samplesheet contains row with ${row.size()} column(s). Expects 5."
                    System.exit(1)
                }
            }
}

// Validate input parameters
WorkflowMag.initialise(params, log, hybrid)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.phix_reference, params.host_fasta, params.kraken2_db, params.gtdb, params.lambda_reference, params.busco_reference ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
                                     CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
                          IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
//def modules = params.modules.clone()
//def multiqc_options   = modules['multiqc']
//multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
//-------------------- MODULE: Local to the pipeline nf-core/mag version 2.1.1 -------------------- 
//
include { CUSTOM_DUMPSOFTWAREVERSIONS as CUSTOM_DUMPSOFTWAREVERSIONS_CUTADAPT } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS as CUSTOM_DUMPSOFTWAREVERSIONS_TRIMMOMATIC } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_HOST_REMOVAL_BUILD                 } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_PHIX_REMOVAL_BUILD_CUTADAPT        } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_PHIX_REMOVAL_BUILD_TRIMMOMATIC     } from '../modules/local/bowtie2_removal_build'
include { KRAKEN2_DB_PREPARATION as KRAKEN2_DB_PREPARATION_CUTADAPT           } from '../modules/local/kraken2_db_preparation'
include { KRONA_DB as KRONA_DB_CUTADAPT                                       } from '../modules/local/krona_db'
include { KRAKEN2_DB_PREPARATION as KRAKEN2_DB_PREPARATION_TRIMMOMATIC        } from '../modules/local/kraken2_db_preparation'
include { KRONA_DB as KRONA_DB_TRIMMOMATIC                                    } from '../modules/local/krona_db'
include { KRAKEN2_DB_PREPARATION                                              } from '../modules/local/kraken2_db_preparation'
include { KRONA_DB                                                            } from '../modules/local/krona_db'
include { POOL_SINGLE_READS                                                   } from '../modules/local/pool_single_reads'
include { POOL_SINGLE_READS as POOL_SINGLE_READS_CUTADAPT                     } from '../modules/local/pool_single_reads'
include { POOL_SINGLE_READS as POOL_SINGLE_READS_TRIMMOMATIC                  } from '../modules/local/pool_single_reads'
include { POOL_PAIRED_READS as POOL_PAIRED_READS_CUTADAPT                     } from '../modules/local/pool_paired_reads'
include { POOL_PAIRED_READS as POOL_PAIRED_READS_TRIMMOMATIC                  } from '../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS                                } from '../modules/local/pool_single_reads'
include { MULTIQC as MULTIQC_CUTADAPT                                         } from '../modules/local/multiqc'                     //addParams( options: multiqc_options                       )
include { MULTIQC as MULTIQC_TRIMMOMATIC                                      } from '../modules/local/multiqc'                     //addParams( options: multiqc_options                       )

//
//-------------------- MODULE: Local to the pipeline edited (or) added by Zifo --------------------
//
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN_CUTADAPT        } from '../modules/local/bowtie2_removal_align'       //addParams( options: modules['bowtie2_host_removal_align_cutadapt'] )
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN_TRIMMOMATIC     } from '../modules/local/bowtie2_removal_align'       //addParams( options: modules['bowtie2_host_removal_align_trimmomatic'] )
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_PHIX_REMOVAL_ALIGN_CUTADAPT        } from '../modules/local/bowtie2_removal_align'       //addParams( options: modules['bowtie2_phix_removal_align_cutadapt'] )
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_PHIX_REMOVAL_ALIGN_TRIMMOMATIC     } from '../modules/local/bowtie2_removal_align'       //addParams( options: modules['bowtie2_phix_removal_align_trimmomatic'] )
include { KRAKEN2 as KRAKEN2_CUTADAPT                                         } from '../modules/local/kraken2'                     //addParams( options: modules['kraken2_cutadapt']           )
include { KRAKEN2 as KRAKEN2_TRIMMOMATIC                                      } from '../modules/local/kraken2'                     //addParams( options: modules['kraken2_trimmomatic']        )
include { KRONA as KRONA_CUTADAPT                                             } from '../modules/local/krona'                       //addParams( options: modules['krona_cutadapt']             )
include { KRONA as KRONA_TRIMMOMATIC     	                                  } from '../modules/local/krona'                       //addParams( options: modules['krona_trimmomatic']          )
include { MEGAHIT as MEGAHIT_CUTADAPT                                         } from '../modules/local/megahit'                     //addParams( options: modules['megahit_cutadapt']           )
include { MEGAHIT as MEGAHIT_TRIMMOMATIC                                      } from '../modules/local/megahit'                     //addParams( options: modules['megahit_trimmomatic']        )
include { SPADES as SPADES_CUTADAPT                                           } from '../modules/local/spades'                      //addParams( options: modules['spades_cutadapt']            )
include { SPADES as SPADES_TRIMMOMATIC                                        } from '../modules/local/spades'                      //addParams( options: modules['spades_trimmomatic']         )
include { QUAST_SPADES as QUAST_CUTADAPT_SPADES                               } from '../modules/local/quast_spades'                //addParams( options: modules['quast_cutadapt']             )
include { QUAST_BINS as QUAST_BINS_CUTADAPT_SPADES                            } from '../modules/local/quast_bins'                  //addParams( options: modules['quast_bins_cutadapt_spades'] )
include { QUAST_BINS_SUMMARY as QUAST_BINS_SUMMARY_CUTADAPT_SPADES            } from '../modules/local/quast_bins_summary'          //addParams( options: modules['quast_bins_summary_cutadapt_spades'])
include { QUAST_SPADES as QUAST_TRIMMOMATIC_SPADES                            } from '../modules/local/quast_spades'                //addParams( options: modules['quast_trimmomatic']          )
include { QUAST_BINS as QUAST_BINS_TRIMMOMATIC_SPADES                         } from '../modules/local/quast_bins'                  //addParams( options: modules['quast_bins_trimmomatic_spades']     )
include { QUAST_BINS_SUMMARY as QUAST_BINS_SUMMARY_TRIMMOMATIC_SPADES         } from '../modules/local/quast_bins_summary'          //addParams( options: modules['quast_bins_summary_trimmomatic_spades']         )

include { QUAST_MEGAHIT as QUAST_CUTADAPT_MEGAHIT                             } from '../modules/local/quast_megahit'               //addParams( options: modules['quast_cutadapt']             )
include { QUAST_BINS as QUAST_BINS_CUTADAPT_MEGAHIT                           } from '../modules/local/quast_bins'                //  addParams( options: modules['quast_bins_cutadapt_megahit'])
include { QUAST_BINS_SUMMARY as QUAST_BINS_SUMMARY_CUTADAPT_MEGAHIT           } from '../modules/local/quast_bins_summary'         // addParams( options: modules['quast_bins_summary_cutadapt_megahit'])
include { QUAST_MEGAHIT as QUAST_TRIMMOMATIC_MEGAHIT                          } from '../modules/local/quast_megahit'             //  addParams( options: modules['quast_trimmomatic']          )
include { QUAST_BINS as QUAST_BINS_TRIMMOMATIC_MEGAHIT                        } from '../modules/local/quast_bins'                 // addParams( options: modules['quast_bins_trimmomatic_megahit']     )
include { QUAST_BINS_SUMMARY as QUAST_BINS_SUMMARY_TRIMMOMATIC_MEGAHIT        } from '../modules/local/quast_bins_summary'         // addParams( options: modules['quast_bins_summary_trimmomatic_megahit']         )
include { BIN_SUMMARY as BIN_SUMMARY_CUTADAPT_SPADES                          } from '../modules/local/bin_summary'                // addParams( options: modules['bin_summary_cutadapt_spades']                )
include { BIN_SUMMARY as BIN_SUMMARY_TRIMMOMATIC_SPADES                       } from '../modules/local/bin_summary'               //  addParams( options: modules['bin_summary_trimmomatic_spades']                )
include { BIN_SUMMARY as BIN_SUMMARY_CUTADAPT_MEGAHIT                         } from '../modules/local/bin_summary'                 //addParams( options: modules['bin_summary_cutadapt_megahit']                )
include { BIN_SUMMARY as BIN_SUMMARY_TRIMMOMATIC_MEGAHIT                      } from '../modules/local/bin_summary'                 //addParams( options: modules['bin_summary_trimmomatic_megahit']                )
include { CUTADAPT                                                            } from '../modules/local/cutadapt'                    //addParams( options: modules['cutadapt']                   )
include { TRIMMOMATIC                                                         } from '../modules/local/trimmomatic'                 //addParams( options: modules['trimmomatic']                )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                                                         } from '../subworkflows/local/input_check'
include { METABAT2_BINNING as METABAT2_BINNING_CUTADAPT_MEGAHIT               } from '../subworkflows/local/metabat2_binning'     //addParams( bowtie2_align_options: modules['bowtie2_assembly_align'], metabat2_options: modules['metabat2_cutadapt'], mag_depths_options: modules['mag_depths_cutadapt'], mag_depths_plot_options: modules['mag_depths_plot'], mag_depths_summary_options: modules['mag_depths_summary_cutadapt_megahit']         )
include { METABAT2_BINNING as METABAT2_BINNING_CUTADAPT_SPADES                } from '../subworkflows/local/metabat2_binning'      //addParams( bowtie2_align_options: modules['bowtie2_assembly_align'], metabat2_options: modules['metabat2_cutadapt'], mag_depths_options: modules['mag_depths_cutadapt'], mag_depths_plot_options: modules['mag_depths_plot'], mag_depths_summary_options: modules['mag_depths_summary_cutadapt_spades']          )
include { METABAT2_BINNING as METABAT2_BINNING_TRIMMOMATIC_MEGAHIT            } from '../subworkflows/local/metabat2_binning'      //addParams( bowtie2_align_options: modules['bowtie2_assembly_align'], metabat2_options: modules['metabat2_trimmomatic'], mag_depths_options: modules['mag_depths_trimmomatic'], mag_depths_plot_options: modules['mag_depths_plot'], mag_depths_summary_options: modules['mag_depths_summary_trimmomatic_megahit'])
include { METABAT2_BINNING as METABAT2_BINNING_TRIMMOMATIC_SPADES             } from '../subworkflows/local/metabat2_binning'      //addParams( bowtie2_align_options: modules['bowtie2_assembly_align'], metabat2_options: modules['metabat2_trimmomatic'], mag_depths_options: modules['mag_depths_trimmomatic'], mag_depths_plot_options: modules['mag_depths_plot'], mag_depths_summary_options: modules['mag_depths_summary_trimmomatic_spades'] )
include { BUSCO_QC as BUSCO_QC_CUTADAPT_SPADES                                } from '../subworkflows/local/busco_qc'              //addParams( busco_db_options: modules['busco_db_preparation'], busco_options: modules['busco_cutadapt_spades'], busco_save_download_options: modules['busco_save_download'], busco_plot_options: modules['busco_plot_cutadapt'], busco_summary_options: modules['busco_summary_cutadapt_spades']                  )
include { BUSCO_QC as BUSCO_QC_TRIMMOMATIC_SPADES                             } from '../subworkflows/local/busco_qc'              // addParams( busco_db_options: modules['busco_db_preparation'], busco_options: modules['busco_trimmomatic_spades'], busco_save_download_options: modules['busco_save_download'], busco_plot_options: modules['busco_plot_trimmomatic'], busco_summary_options: modules['busco_summary_trimmomatic_spades']         )
include { BUSCO_QC as BUSCO_QC_CUTADAPT_MEGAHIT                               } from '../subworkflows/local/busco_qc'              // addParams( busco_db_options: modules['busco_db_preparation'], busco_options: modules['busco_cutadapt_megahit'], busco_save_download_options: modules['busco_save_download'], busco_plot_options: modules['busco_plot_cutadapt'], busco_summary_options: modules['busco_summary_cutadapt_megahit']                )
include { BUSCO_QC as BUSCO_QC_TRIMMOMATIC_MEGAHIT                            } from '../subworkflows/local/busco_qc'              // addParams( busco_db_options: modules['busco_db_preparation'], busco_options: modules['busco_trimmomatic_megahit'], busco_save_download_options: modules['busco_save_download'], busco_plot_options: modules['busco_plot_trimmomatic'], busco_summary_options: modules['busco_summary_trimmomatic_megahit']       )
include { GTDBTK                                        } from '../subworkflows/local/gtdbtk'                //addParams( gtdbtk_classify_options: modules['gtdbtk_classify'], gtdbtk_summary_options: modules['gtdbtk_summary'])
include { SOURMASH as SOURMASH_CUTADAPT_SPADES                                } from '../subworkflows/local/sourmash'              // addParams( options: modules['sourmash_cutadapt_spades']                    )
include { SOURMASH as SOURMASH_CUTADAPT_MEGAHIT                               } from '../subworkflows/local/sourmash'               //addParams( options: modules['sourmash_cutadapt_megahit']                   )
include { SOURMASH as SOURMASH_TRIMMOMATIC_SPADES                             } from '../subworkflows/local/sourmash'               //addParams( options: modules['sourmash_trimmomatic_spades']                 )
include { SOURMASH as SOURMASH_TRIMMOMATIC_MEGAHIT                            } from '../subworkflows/local/sourmash'             //  addParams( options: modules['sourmash_trimmomatic_megahit']                )
include { GTDBTK  as GTDBTK_CUTADAPT_MEGAHIT                                  } from '../subworkflows/local/gtdbtk'               //  addParams( gtdbtk_classify_options: modules['gtdbtk_classify'], gtdbtk_summary_options: modules['gtdbtk_summary_cutadapt_megahit']      )
include { GTDBTK  as GTDBTK_CUTADAPT_SPADES                                   } from '../subworkflows/local/gtdbtk'                // addParams( gtdbtk_classify_options: modules['gtdbtk_classify'], gtdbtk_summary_options: modules['gtdbtk_summary_cutadapt_spades']       )
include { GTDBTK  as GTDBTK_TRIMMOMATIC_MEGAHIT                               } from '../subworkflows/local/gtdbtk'                // addParams( gtdbtk_classify_options: modules['gtdbtk_classify'], gtdbtk_summary_options: modules['gtdbtk_summary_trimmomatic_megahit']   )
include { GTDBTK  as GTDBTK_TRIMMOMATIC_SPADES                                } from '../subworkflows/local/gtdbtk'                // addParams( gtdbtk_classify_options: modules['gtdbtk_classify'], gtdbtk_summary_options: modules['gtdbtk_summary_trimmomatic_spades']    )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
//---------------- MODULE: Installed directly from nf-core/modules ---------------------
//
include { FASTQC as FASTQC_RAW                                                } from '../modules/nf-core/modules/fastqc/main'       //addParams( options: modules['fastqc_raw']            )
include { FASTQC as FASTQC_TRIMMED_CUTADAPT                                   } from '../modules/nf-core/modules/fastqc/main'       //addParams( options: modules['fastqc_trimmed_cutadapt'])
include { FASTQC as FASTQC_TRIMMED_TRIMMOMATIC                                } from '../modules/nf-core/modules/fastqc/main'       //addParams( options: modules['fastqc_trimmed_trimmomatic'])
include { PRODIGAL as PRODIGAL_SPADES_CUTADAPT                                } from '../modules/nf-core/modules/prodigal/main'     //addParams( options: modules['prodigal']              )
include { PRODIGAL as PRODIGAL_MEGAHIT_CUTADAPT                               } from '../modules/nf-core/modules/prodigal/main'     //addParams( options: modules['prodigal']              )
include { PRODIGAL as PRODIGAL_SPADES_TRIMMOMATIC                             } from '../modules/nf-core/modules/prodigal/main'     //addParams( options: modules['prodigal']              )
include { PRODIGAL as PRODIGAL_MEGAHIT_TRIMMOMATIC                            } from '../modules/nf-core/modules/prodigal/main'     //addParams( options: modules['prodigal']              )


//
//---------------- MODULE: Installed directly from nf-core/modules edited (or) added by Zifo ----------------
//
include { PROKKA as PROKKA_BINS_CUTADAPT_SPADES                               } from '../modules/nf-core/modules/prokka/main'       //addParams( options: modules['prokka']                )
include { PROKKA as PROKKA_BINS_CUTADAPT_MEGAHIT                              } from '../modules/nf-core/modules/prokka/main'       //addParams( options: modules['prokka']                )
include { PROKKA as PROKKA_BINS_TRIMMOMATIC_SPADES                            } from '../modules/nf-core/modules/prokka/main'       //addParams( options: modules['prokka']                )
include { PROKKA as PROKKA_BINS_TRIMMOMATIC_MEGAHIT                           } from '../modules/nf-core/modules/prokka/main'       //addParams( options: modules['prokka']                )
include { PROKKA as PROKKA_ASSEMBLIES_CUTADAPT_SPADES                         } from '../modules/nf-core/modules/prokka/main'       //addParams( options: modules['prokka']                )
include { PROKKA as PROKKA_ASSEMBLIES_CUTADAPT_MEGAHIT                        } from '../modules/nf-core/modules/prokka/main'       //addParams( options: modules['prokka']                )
include { PROKKA as PROKKA_ASSEMBLIES_TRIMMOMATIC_SPADES                      } from '../modules/nf-core/modules/prokka/main'       //addParams( options: modules['prokka']                )
include { PROKKA as PROKKA_ASSEMBLIES_TRIMMOMATIC_MEGAHIT                     } from '../modules/nf-core/modules/prokka/main'       //addParams( options: modules['prokka']                )

////////////////////////////////////////////////////
/* --  Create channel for reference databases  -- */
////////////////////////////////////////////////////

if ( params.host_genome ) {
    host_fasta = params.genomes[params.host_genome].fasta ?: false
    Channel
        .value(file( "${host_fasta}" ))
        .set { ch_host_fasta }

    host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
    Channel
        .value(file( "${host_bowtie2index}/*" ))
        .set { ch_host_bowtie2index }
} else if ( params.host_fasta ) {
    Channel
        .value(file( "${params.host_fasta}" ))
        .set { ch_host_fasta }
} else {
    ch_host_fasta = Channel.empty()
}

if(params.busco_reference){
    Channel
        .value(file( "${params.busco_reference}" ))
        .set { ch_busco_db_file }
} else {
    ch_busco_db_file = Channel.empty()
}
if (params.busco_download_path) {
    Channel
        .value(file( "${params.busco_download_path}" ))
        .set { ch_busco_download_folder }
} else {
    ch_busco_download_folder = Channel.empty()
}

if(params.kraken2_db){
    Channel
        .value(file( "${params.kraken2_db}" ))
        .set { ch_kraken2_db_file }
} else {
    ch_kraken2_db_file = Channel.empty()
}

if(!params.keep_phix) {
    Channel
        .value(file( "${params.phix_reference}" ))
        .set { ch_phix_db_file }
}

if (!params.keep_lambda) {
    Channel
        .value(file( "${params.lambda_reference}" ))
        .set { ch_nanolyse_db }
}

gtdb = params.skip_busco ? false : params.gtdb
if (params.gtdbtk && gtdb) {
    Channel
        .value(file( "${gtdb}" ))
        .set { ch_gtdb }
} else {
    ch_gtdb = Channel.empty()
}


//-----------------------------------Added by Zifo--------------------------------------

if(params.sourmash_db){
    Channel
        .value(file( "${params.sourmash_db}" ))
        .set { ch_sourmash_db }
} else {
    ch_sourmash_db = Channel.empty()
}


/*
========================================================================================
                                  RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report    = []
def busco_failed_bins = [:]

workflow METABP {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ()
    ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
    ch_raw_long_reads  = INPUT_CHECK.out.raw_long_reads

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    FASTQC_RAW (
        ch_raw_short_reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_RAW.out.versions.first().ifEmpty(null))

    /*
    ---------------------Cutadapt added for Quality Control by Zifo-------------------
    */
    if (params.cutadapt) {
    
    CUTADAPT (
        ch_raw_short_reads
    )
    CUTADAPT.out.reads
            .map { meta, reads ->
                def meta_new = meta.clone()
                meta_new.trimmer  = "Cutadapt"
                [ meta_new, reads ]
            }
            .set { ch_short_reads_cutadapt }
    
    ch_software_versions = ch_software_versions.mix(CUTADAPT.out.versions.first().ifEmpty(null))
    
    ch_bowtie2_removal_host_multiqc = Channel.empty()
    if (params.host_fasta || params.host_genome) {
        BOWTIE2_HOST_REMOVAL_ALIGN_CUTADAPT (
            ch_short_reads_cutadapt,
            ch_host_bowtie2index
        )

        ch_short_reads_cutadapt = BOWTIE2_HOST_REMOVAL_ALIGN_CUTADAPT.out.reads
        ch_bowtie2_removal_host_multiqc = BOWTIE2_HOST_REMOVAL_ALIGN_CUTADAPT.out.log 
        ch_software_versions = ch_software_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN_CUTADAPT.out.versions.first().ifEmpty(null))
    }
    
    if(!params.keep_phix) {
        BOWTIE2_PHIX_REMOVAL_BUILD_CUTADAPT (
            ch_phix_db_file
        )
        BOWTIE2_PHIX_REMOVAL_ALIGN_CUTADAPT (
            ch_short_reads_cutadapt,
            BOWTIE2_PHIX_REMOVAL_BUILD_CUTADAPT.out.index
        )
	
        ch_short_reads_cutadapt = BOWTIE2_PHIX_REMOVAL_ALIGN_CUTADAPT.out.reads 
        ch_short_reads = ch_short_reads_cutadapt
    }

    FASTQC_TRIMMED_CUTADAPT (
        ch_short_reads_cutadapt
    )

    ch_fastqc_trimmed = FASTQC_TRIMMED_CUTADAPT.out.zip

    if (params.kraken2) {
        KRAKEN2_DB_PREPARATION_CUTADAPT (
        ch_kraken2_db_file
        )
        KRAKEN2_CUTADAPT (
            ch_short_reads_cutadapt,
            KRAKEN2_DB_PREPARATION_CUTADAPT.out.db
        )

        ch_software_versions = ch_software_versions.mix(KRAKEN2_CUTADAPT.out.versions.first().ifEmpty(null))

        if ( params.kraken2_db && !params.skip_krona){
            KRONA_DB_CUTADAPT ()
            KRAKEN2_CUTADAPT.out.results_for_krona
                . map { classifier, meta, report ->
                    def meta_new = meta.clone()
                    meta_new.classifier  = classifier
                    [ meta_new, report ]
                }
                .set { ch_tax_classifications_cutadapt }
        
            KRONA_CUTADAPT (
                ch_tax_classifications_cutadapt,
                KRONA_DB_CUTADAPT.out.db.collect()
            )
            ch_software_versions = ch_software_versions.mix(KRONA_CUTADAPT.out.versions.first().ifEmpty(null))
        }
    }

    /*
    ================================================================================
                                      Assembly
    ================================================================================
    */

    // Co-assembly: prepare grouping for MEGAHIT and for pooling for SPAdes
    if (params.coassemble_group) {
        // short reads
        // group and set group as new id
        ch_short_reads_cutadapt
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                    def meta = [:]
                    meta.id          = "group-$group"
                    meta.group       = group
                    meta.single_end  = params.single_end
                    if (!params.single_end) [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
                    else [ meta, reads.collect { it }, [] ]
            }
            .set { ch_short_reads_grouped_cutadapt }

    } else {
        ch_short_reads_cutadapt
            .map { meta, reads ->
                    if (!params.single_end){ [ meta, [reads[0]], [reads[1]] ] }
                    else [ meta, [reads], [] ] }
            .set { ch_short_reads_grouped_cutadapt }
    }


    /*
    -------------------------------------------------------------------------------------------------
                                            Assembly of reads
    -------------------------------------------------------------------------------------------------
    */

    ch_assemblies_cutadapt = Channel.empty()
    if (!params.skip_megahit) {
        MEGAHIT_CUTADAPT ( ch_short_reads_grouped_cutadapt )
        MEGAHIT_CUTADAPT.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "MEGAHIT"
                [ meta_new, assembly ]
            }
            .set { ch_megahit_assemblies_cutadapt }
        ch_assemblies_cutadapt = ch_assemblies_cutadapt.mix(ch_megahit_assemblies_cutadapt)
        ch_software_versions = ch_software_versions.mix(MEGAHIT_CUTADAPT.out.versions.first().ifEmpty(null))
    }

    // Co-assembly: pool reads for SPAdes
    if (params.coassemble_group) {
        // short reads
        if (!params.single_end && (!params.skip_spades || !params.skip_spadeshybrid)){
            if (params.single_end){
		        /*POOL_SINGLE_READS ( ch_short_reads_grouped )
                ch_short_reads_spades = POOL_SINGLE_READS.out.reads*/
                POOL_SINGLE_READS_CUTADAPT ( ch_short_reads_grouped_cutadapt )
                ch_short_reads_spades_cutadapt = POOL_SINGLE_READS_CUTADAPT.out.reads
            } else {
                POOL_PAIRED_READS_CUTADAPT ( ch_short_reads_grouped_cutadapt )
                ch_short_reads_spades_cutadapt = POOL_PAIRED_READS_CUTADAPT.out.reads
            }
        }
    } else {
        ch_short_reads_spades_cutadapt = ch_short_reads_cutadapt
        }
	
	    //ch_short_reads_spades_cutadapt = ch_short_reads_cutadapt

    if (!params.single_end && !params.skip_spades){
        SPADES_CUTADAPT ( ch_short_reads_spades_cutadapt )
        SPADES_CUTADAPT.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdes"
                [ meta_new, assembly ]
            }
            .set { ch_spades_assemblies_cutadapt }
	    SPADES_CUTADAPT.out.contig
            .map { meta, contig ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdes"
                [ meta_new, contig ]
            }
            .set { ch_spades_contigs_cutadapt }
        ch_assemblies_cutadapt = Channel.empty()
	    ch_assemblies_cutadapt = ch_assemblies_cutadapt.mix(ch_spades_assemblies_cutadapt)
        ch_software_versions = ch_software_versions.mix(SPADES_CUTADAPT.out.versions.first().ifEmpty(null))
    }

    /*
    -------------------------------------------------------------------------------------------------
                                                Assembly QC
    -------------------------------------------------------------------------------------------------
    */
	
    ch_quast_multiqc = Channel.empty()
    if (!params.skip_quast){
        if (!params.single_end && !params.skip_spades) {
            QUAST_CUTADAPT_SPADES ( ch_spades_assemblies_cutadapt, ch_spades_contigs_cutadapt )
        	ch_quast_multiqc = ch_quast_multiqc.mix(QUAST_CUTADAPT_SPADES.out.qc)
            ch_software_versions = ch_software_versions.mix(QUAST_CUTADAPT_SPADES.out.versions.first().ifEmpty(null))
        }
        if (!params.skip_megahit){
	        QUAST_CUTADAPT_MEGAHIT ( ch_megahit_assemblies_cutadapt )
            ch_quast_multiqc = ch_quast_multiqc.mix(QUAST_CUTADAPT_MEGAHIT.out.qc)
            ch_software_versions = ch_software_versions.mix(QUAST_CUTADAPT_MEGAHIT.out.versions.first().ifEmpty(null))
        }
    }
    
    /*
    ================================================================================================
                                            Predict proteins
    ================================================================================================
    */

    if (!params.skip_prodigal){
        if (!params.single_end && !params.skip_spades) {
            PRODIGAL_SPADES_CUTADAPT (
            ch_spades_assemblies_cutadapt,
			'gff')
           // modules['prodigal']['output_format']
       
        ch_software_versions = ch_software_versions.mix(PRODIGAL_SPADES_CUTADAPT.out.versions.first().ifEmpty(null))
            }
        if (!params.skip_megahit){
	        PRODIGAL_MEGAHIT_CUTADAPT (
            ch_megahit_assemblies_cutadapt,
			'gff')
            //modules['prodigal']['output_format']
       
        ch_software_versions = ch_software_versions.mix(PRODIGAL_MEGAHIT_CUTADAPT.out.versions.first().ifEmpty(null))
        }
    }
    
    /*
    ================================================================================================
                                                Binning
    ================================================================================================
    */

    ch_bowtie2_assembly_multiqc = Channel.empty()
    ch_busco_summary            = Channel.empty()
    ch_busco_multiqc            = Channel.empty()
    if (!params.skip_binning) {
        if (!params.single_end && !params.skip_spades) {
            METABAT2_BINNING_CUTADAPT_SPADES (
            ch_spades_assemblies_cutadapt,
            ch_short_reads_cutadapt
        )
        ch_bowtie2_assembly_multiqc = ch_bowtie2_assembly_multiqc.mix(METABAT2_BINNING_CUTADAPT_SPADES.out.bowtie2_assembly_multiqc)
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING_CUTADAPT_SPADES.out.bowtie2_versions.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING_CUTADAPT_SPADES.out.metabat2_versions.first().ifEmpty(null))
        }
        if (!params.skip_megahit) {
            METABAT2_BINNING_CUTADAPT_MEGAHIT (
            ch_megahit_assemblies_cutadapt,
            ch_short_reads_cutadapt
        )
        ch_bowtie2_assembly_multiqc = ch_bowtie2_assembly_multiqc.mix(METABAT2_BINNING_CUTADAPT_MEGAHIT.out.bowtie2_assembly_multiqc)
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING_CUTADAPT_MEGAHIT.out.bowtie2_versions.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING_CUTADAPT_MEGAHIT.out.metabat2_versions.first().ifEmpty(null))
        }

        /*
         * BUSCO subworkflow: Quantitative measures for the assessment of genome assembly
         */
        if (!params.skip_busco) {
            if (!params.single_end && !params.skip_spades) {
                BUSCO_QC_CUTADAPT_SPADES (
                    ch_busco_db_file,
                    ch_busco_download_folder,
                    METABAT2_BINNING_CUTADAPT_SPADES.out.bins.transpose()
                )
                ch_busco_summary_cutadapt_spades = BUSCO_QC_CUTADAPT_SPADES.out.summary
                ch_busco_summary = ch_busco_summary.mix(ch_busco_summary_cutadapt_spades)
                ch_busco_multiqc = ch_busco_multiqc.mix(BUSCO_QC_CUTADAPT_SPADES.out.multiqc)
                ch_software_versions = ch_software_versions.mix(BUSCO_QC_CUTADAPT_SPADES.out.versions.first().ifEmpty(null))
            
                // process information if BUSCO analysis failed for individual bins due to no matching genes
                BUSCO_QC_CUTADAPT_SPADES.out
                    .failed_bin
                    .splitCsv(sep: '\t')
                    .map { bin, error -> if (!bin.contains(".unbinned.")) busco_failed_bins[bin] = error }
            }
            if (!params.skip_megahit) {
                BUSCO_QC_CUTADAPT_MEGAHIT (
                    ch_busco_db_file,
                    ch_busco_download_folder,
                    METABAT2_BINNING_CUTADAPT_MEGAHIT.out.bins.transpose()
                )
                ch_busco_summary_cutadapt_megahit = BUSCO_QC_CUTADAPT_MEGAHIT.out.summary
                ch_busco_summary = ch_busco_summary.mix(ch_busco_summary_cutadapt_megahit)
                ch_busco_multiqc = ch_busco_multiqc.mix(BUSCO_QC_CUTADAPT_MEGAHIT.out.multiqc)
                ch_software_versions = ch_software_versions.mix(BUSCO_QC_CUTADAPT_MEGAHIT.out.versions.first().ifEmpty(null))
            
                // process information if BUSCO analysis failed for individual bins due to no matching genes
                BUSCO_QC_CUTADAPT_MEGAHIT.out
                    .failed_bin
                    .splitCsv(sep: '\t')
                    .map { bin, error -> if (!bin.contains(".unbinned.")) busco_failed_bins[bin] = error }
            }
        }

        if (!params.skip_quast){
            if (!params.single_end && !params.skip_spades) {
                QUAST_BINS_CUTADAPT_SPADES ( METABAT2_BINNING_CUTADAPT_SPADES.out.bins )
                ch_software_versions = ch_software_versions.mix(QUAST_BINS_CUTADAPT_SPADES.out.versions.first().ifEmpty(null))

                QUAST_BINS_SUMMARY_CUTADAPT_SPADES ( QUAST_BINS_CUTADAPT_SPADES.out.quast_bin_summaries.collect() )
                ch_quast_bins_summary_cutadapt_spades = QUAST_BINS_SUMMARY_CUTADAPT_SPADES.out.summary
            }
            if (!params.skip_megahit) {
                QUAST_BINS_CUTADAPT_MEGAHIT ( METABAT2_BINNING_CUTADAPT_MEGAHIT.out.bins )
                ch_software_versions = ch_software_versions.mix(QUAST_BINS_CUTADAPT_MEGAHIT.out.versions.first().ifEmpty(null))

                QUAST_BINS_SUMMARY_CUTADAPT_MEGAHIT ( QUAST_BINS_CUTADAPT_MEGAHIT.out.quast_bin_summaries.collect() )
                ch_quast_bins_summary_cutadapt_megahit = QUAST_BINS_SUMMARY_CUTADAPT_MEGAHIT.out.summary
            }
        }

        /*
        -------------------------------------------------------------------------------------------------
                                Sourmash added for Taxonomic Classification
        -------------------------------------------------------------------------------------------------
        */

        /*
         * Sourmash subworkflow: k-mer based taxonomic exploration and classification routines for genome and metagenome analysis
         */
        if (params.sourmash) {
            /* if (!params.single_end && !params.skip_spades) {
                SOURMASH_SIGNATURE_CUTADAPT_SPADES ( METABAT2_BINNING_CUTADAPT_SPADES.out.bins )
                
                SOURMASH_SUMMARIZE_CUTADAPT_SPADES ( 
                    ch_sourmash_db,
                    SOURMASH_SIGNATURE_CUTADAPT_SPADES.out.signatures
                )
                ch_software_versions = ch_software_versions.mix(SOURMASH_SIGNATURE_CUTADAPT_SPADES.out.versions.first().ifEmpty(null))
            }
            if (!params.skip_megahit) {
                SOURMASH_SIGNATURE_CUTADAPT_MEGAHIT ( METABAT2_BINNING_CUTADAPT_MEGAHIT.out.bins )
                
                SOURMASH_SUMMARIZE_CUTADAPT_MEGAHIT ( 
                    ch_sourmash_db,
                    SOURMASH_SIGNATURE_CUTADAPT_MEGAHIT.out.signatures
                )
                ch_software_versions = ch_software_versions.mix(SOURMASH_SIGNATURE_CUTADAPT_MEGAHIT.out.versions.first().ifEmpty(null))
            } */

            if (!params.single_end && !params.skip_spades) {
                SOURMASH_CUTADAPT_SPADES (ch_sourmash_db, METABAT2_BINNING_CUTADAPT_SPADES.out.bins)
                ch_software_versions = ch_software_versions.mix(SOURMASH_CUTADAPT_SPADES.out.versions.first().ifEmpty(null))
            }
            if (!params.skip_megahit) {
                SOURMASH_CUTADAPT_MEGAHIT ( ch_sourmash_db, METABAT2_BINNING_CUTADAPT_MEGAHIT.out.bins )
                ch_software_versions = ch_software_versions.mix(SOURMASH_CUTADAPT_MEGAHIT.out.versions.first().ifEmpty(null))
            }
        }

        /*
         * GTDB-tk: taxonomic classifications using GTDB reference
         */
        if (params.gtdbtk && gtdb) {
            ch_gtdbtk_summary_cutadapt_spades = Channel.empty()
            ch_gtdbtk_summary_cutadapt_megahit = Channel.empty()
            if (!params.single_end && !params.skip_spades) {
                GTDBTK_CUTADAPT_SPADES ( 
                    METABAT2_BINNING_CUTADAPT_SPADES.out.bins,
                    ch_busco_summary_cutadapt_spades,
                    ch_gtdb
                )
                ch_software_versions = ch_software_versions.mix(GTDBTK_CUTADAPT_SPADES.out.versions.first().ifEmpty(null))
                ch_gtdbtk_summary_cutadapt_spades = GTDBTK_CUTADAPT_SPADES.out.summary
                if (!params.skip_busco || !params.skip_quast || gtdb) {
                    BIN_SUMMARY_CUTADAPT_SPADES (
                        METABAT2_BINNING_CUTADAPT_SPADES.out.depths_summary,
                        ch_busco_summary_cutadapt_spades.ifEmpty([]),
                        ch_quast_bins_summary_cutadapt_spades.ifEmpty([]),
                        ch_gtdbtk_summary_cutadapt_spades.ifEmpty([])
                    )
                    ch_bin_depths_summary = METABAT2_BINNING_CUTADAPT_SPADES.out.depths_summary
                    ch_busco_summary = ch_busco_summary_cutadapt_spades
                    ch_quast_bins_summary = ch_quast_bins_summary_cutadapt_spades
                }
            }
            if (!params.skip_megahit) {
                GTDBTK_CUTADAPT_MEGAHIT ( 
                    METABAT2_BINNING_CUTADAPT_MEGAHIT.out.bins,
                    ch_busco_summary_cutadapt_megahit,
                    ch_gtdb
                )
                ch_software_versions = ch_software_versions.mix(GTDBTK_CUTADAPT_MEGAHIT.out.versions.first().ifEmpty(null))
                ch_gtdbtk_summary_cutadapt_megahit = GTDBTK_CUTADAPT_MEGAHIT.out.summary
                if (!params.skip_busco || !params.skip_quast || gtdb) {
                    BIN_SUMMARY_CUTADAPT_MEGAHIT (
                        METABAT2_BINNING_CUTADAPT_MEGAHIT.out.depths_summary,
                        ch_busco_summary_cutadapt_megahit.ifEmpty([]),
                        ch_quast_bins_summary_cutadapt_megahit.ifEmpty([]),
                        ch_gtdbtk_summary_cutadapt_megahit.ifEmpty([])
                    )
                    ch_bin_depths_summary = METABAT2_BINNING_CUTADAPT_MEGAHIT.out.depths_summary
                    ch_busco_summary = ch_busco_summary_cutadapt_megahit
                    ch_quast_bins_summary = ch_quast_bins_summary_cutadapt_megahit
                }
            }
        }

        /*
         * Prokka: Genome annotation
         */
        if (!params.single_end && !params.skip_spades) {
            METABAT2_BINNING_CUTADAPT_SPADES.out.bins
                .transpose()
                    .map { meta, bin ->
                        def meta_new = meta.clone()
                        meta_new.id  = bin.getBaseName()
                        [ meta_new, bin ]
                    }
                    .set { ch_bins_for_prokka_cutadapt_spades }

            if (!params.skip_prokka){
                PROKKA_BINS_CUTADAPT_SPADES (
                    ch_bins_for_prokka_cutadapt_spades,
                    [],
                    []
                )
                PROKKA_ASSEMBLIES_CUTADAPT_SPADES (
                    ch_spades_assemblies_cutadapt,
                    [],
                    []
                )
                ch_software_versions = ch_software_versions.mix(PROKKA_BINS_CUTADAPT_SPADES.out.versions.first().ifEmpty(null))
            }
        }
        if (!params.skip_megahit) {
            METABAT2_BINNING_CUTADAPT_MEGAHIT.out.bins
                .transpose()
                    .map { meta, bin ->
                        def meta_new = meta.clone()
                        meta_new.id  = bin.getBaseName()
                        [ meta_new, bin ]
                    }
                    .set { ch_bins_for_prokka_cutadapt_megahit }

            if (!params.skip_prokka){
                PROKKA_BINS_CUTADAPT_MEGAHIT (
                    ch_bins_for_prokka_cutadapt_megahit,
                    [],
                    []
                )
                PROKKA_ASSEMBLIES_CUTADAPT_MEGAHIT (
                    ch_megahit_assemblies_cutadapt,
                    [],
                    []
                )
                ch_software_versions = ch_software_versions.mix(PROKKA_BINS_CUTADAPT_MEGAHIT.out.versions.first().ifEmpty(null))
            }
        }
        
    }


    //
    // MODULE: Pipeline reporting
    //
    /* ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS_CUTADAPT (
        ch_software_versions.map { it }.collect()
    ) */

    CUSTOM_DUMPSOFTWAREVERSIONS_CUTADAPT (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMag.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS_CUTADAPT.out.mqc_yml.collect())

    MULTIQC_CUTADAPT (
        ch_multiqc_files.collect(),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIMMED_CUTADAPT.out.zip.collect{it[1]}.ifEmpty([]),
        ch_bowtie2_removal_host_multiqc.collect{it[1]}.ifEmpty([]),
        ch_quast_multiqc.collect().ifEmpty([]),
        ch_bowtie2_assembly_multiqc.collect().ifEmpty([]),
        ch_busco_multiqc.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC_CUTADAPT.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC_CUTADAPT.out.versions.ifEmpty(null))
    
    
    }



    /*
    ---------------------Trimmomatic added for Quality Control by Zifo-------------------
    */
    if (params.trimmomatic) {

    TRIMMOMATIC (
        ch_raw_short_reads
    )
    TRIMMOMATIC.out.reads
            .map { meta, reads ->
                def meta_new = meta.clone()
                meta_new.trimmer  = "Trimmomatic"
                [ meta_new, reads ]
            }
            .set { ch_short_reads_trimmomatic }
    ch_software_versions = ch_software_versions.mix(TRIMMOMATIC.out.versions.first().ifEmpty(null))

    ch_bowtie2_removal_host_multiqc = Channel.empty()
    if (params.host_fasta || params.host_genome) {
        BOWTIE2_HOST_REMOVAL_ALIGN_TRIMMOMATIC (
            ch_short_reads_trimmomatic,
            ch_host_bowtie2index
        )

        ch_short_reads_trimmomatic = BOWTIE2_HOST_REMOVAL_ALIGN_TRIMMOMATIC.out.reads
        ch_bowtie2_removal_host_multiqc = BOWTIE2_HOST_REMOVAL_ALIGN_TRIMMOMATIC.out.log 
        ch_software_versions = ch_software_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN_TRIMMOMATIC.out.versions.first().ifEmpty(null))
    }

    if(!params.keep_phix) {
        BOWTIE2_PHIX_REMOVAL_BUILD_TRIMMOMATIC (
            ch_phix_db_file
        )
        BOWTIE2_PHIX_REMOVAL_ALIGN_TRIMMOMATIC (
            ch_short_reads_trimmomatic,
            BOWTIE2_PHIX_REMOVAL_BUILD_TRIMMOMATIC.out.index
        )

        ch_short_reads_trimmomatic = BOWTIE2_PHIX_REMOVAL_ALIGN_TRIMMOMATIC.out.reads 
        ch_short_reads = ch_short_reads_trimmomatic//.mix(ch_short_reads_trimmomatic)
    }

    FASTQC_TRIMMED_TRIMMOMATIC (
        ch_short_reads_trimmomatic
    )

    ch_fastqc_trimmed = FASTQC_TRIMMED_TRIMMOMATIC.out.zip

    if (params.kraken2) {
        KRAKEN2_DB_PREPARATION (
        ch_kraken2_db_file
        )
        KRAKEN2_TRIMMOMATIC (
            ch_short_reads_trimmomatic,
            KRAKEN2_DB_PREPARATION.out.db
        )

        ch_software_versions = ch_software_versions.mix(KRAKEN2_TRIMMOMATIC.out.versions.first().ifEmpty(null))

        if ( params.kraken2_db && !params.skip_krona){
            KRONA_DB_TRIMMOMATIC ()
            KRAKEN2_TRIMMOMATIC.out.results_for_krona
                . map { classifier, meta, report ->
                    def meta_new = meta.clone()
                    meta_new.classifier  = classifier
                    [ meta_new, report ]
                }
                .set { ch_tax_classifications_trimmomatic }

            KRONA_TRIMMOMATIC (
                ch_tax_classifications_trimmomatic,
                KRONA_DB_TRIMMOMATIC.out.db.collect()
            )
            ch_software_versions = ch_software_versions.mix(KRONA_TRIMMOMATIC.out.versions.first().ifEmpty(null))
        }
    }

    /*
    ================================================================================
                                      Assembly
    ================================================================================
    */

    // Co-assembly: prepare grouping for MEGAHIT and for pooling for SPAdes
    if (params.coassemble_group) {
        // short reads
        // group and set group as new id
        ch_short_reads_trimmomatic
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                    def meta = [:]
                    meta.id          = "group-$group"
                    meta.group       = group
                    meta.single_end  = params.single_end
                    if (!params.single_end) [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
                    else [ meta, reads.collect { it }, [] ]
            }
            .set { ch_short_reads_grouped_trimmomatic }

    } else {
        ch_short_reads_trimmomatic
            .map { meta, reads ->
                    if (!params.single_end){ [ meta, [reads[0]], [reads[1]] ] }
                    else [ meta, [reads], [] ] }
            .set { ch_short_reads_grouped_trimmomatic }
    }


    /*
    -------------------------------------------------------------------------------------------------
                                            Assembly of reads
    -------------------------------------------------------------------------------------------------
    */

    ch_assemblies_trimmomatic = Channel.empty()
    if (!params.skip_megahit) {
        MEGAHIT_TRIMMOMATIC ( ch_short_reads_grouped_trimmomatic )
        MEGAHIT_TRIMMOMATIC.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "MEGAHIT"
                [ meta_new, assembly ]
            }
            .set { ch_megahit_assemblies_trimmomatic }
        ch_assemblies_trimmomatic = ch_assemblies_trimmomatic.mix(ch_megahit_assemblies_trimmomatic)
        ch_software_versions = ch_software_versions.mix(MEGAHIT_TRIMMOMATIC.out.versions.first().ifEmpty(null))
    }

    // Co-assembly: pool reads for SPAdes
    if (params.coassemble_group) {
        // short reads
        if (!params.single_end && (!params.skip_spades || !params.skip_spadeshybrid)){
            if (params.single_end){
		/*POOL_SINGLE_READS ( ch_short_reads_grouped )
                ch_short_reads_spades = POOL_SINGLE_READS.out.reads*/
                POOL_SINGLE_READS_TRIMMOMATIC ( ch_short_reads_grouped_trimmomatic )
                ch_short_reads_spades_trimmomatic = POOL_SINGLE_READS_TRIMMOMATIC.out.reads
            } else {
                POOL_PAIRED_READS_TRIMMOMATIC ( ch_short_reads_grouped_trimmomatic )
                ch_short_reads_spades_trimmomatic = POOL_PAIRED_READS_TRIMMOMATIC.out.reads
            }
        }
      } else {
        ch_short_reads_spades_trimmomatic = ch_short_reads_trimmomatic
      }

    if (!params.single_end && !params.skip_spades){
        SPADES_TRIMMOMATIC ( ch_short_reads_spades_trimmomatic )
        SPADES_TRIMMOMATIC.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdes"
                [ meta_new, assembly ]
            }
            .set { ch_spades_assemblies_trimmomatic }
	    SPADES_TRIMMOMATIC.out.contig
            .map { meta, contig ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdes"
                [ meta_new, contig ]
            }
            .set { ch_spades_contigs_trimmomatic }
        ch_assemblies_trimmomatic = Channel.empty()
	    ch_assemblies_trimmomatic = ch_assemblies_trimmomatic.mix(ch_spades_assemblies_trimmomatic)
        ch_software_versions = ch_software_versions.mix(SPADES_TRIMMOMATIC.out.versions.first().ifEmpty(null))
    }

    /*
    -------------------------------------------------------------------------------------------------
                                                Assembly QC
    -------------------------------------------------------------------------------------------------
    */

    ch_quast_multiqc = Channel.empty()
    if (!params.skip_quast){
        if (!params.single_end && !params.skip_spades) {
            QUAST_TRIMMOMATIC_SPADES ( ch_spades_assemblies_trimmomatic, ch_spades_contigs_trimmomatic )
        	ch_quast_multiqc = ch_quast_multiqc.mix(QUAST_TRIMMOMATIC_SPADES.out.qc)
            ch_software_versions = ch_software_versions.mix(QUAST_TRIMMOMATIC_SPADES.out.versions.first().ifEmpty(null))
        }
        if (!params.skip_megahit){
	        QUAST_TRIMMOMATIC_MEGAHIT ( ch_megahit_assemblies_trimmomatic )
            ch_quast_multiqc = ch_quast_multiqc.mix(QUAST_TRIMMOMATIC_MEGAHIT.out.qc)
            ch_software_versions = ch_software_versions.mix(QUAST_TRIMMOMATIC_MEGAHIT.out.versions.first().ifEmpty(null))
        }
    }

    /*
    ================================================================================================
                                            Predict proteins
    ================================================================================================
    */

    if (!params.skip_prodigal){
        if (!params.single_end && !params.skip_spades) {
            PRODIGAL_SPADES_TRIMMOMATIC (
            ch_spades_assemblies_trimmomatic,
            'gff')
        
        ch_software_versions = ch_software_versions.mix(PRODIGAL_SPADES_TRIMMOMATIC.out.versions.first().ifEmpty(null))
            }
        if (!params.skip_megahit){
	        PRODIGAL_MEGAHIT_TRIMMOMATIC (
            ch_megahit_assemblies_trimmomatic,
            'gff')
        
        ch_software_versions = ch_software_versions.mix(PRODIGAL_MEGAHIT_TRIMMOMATIC.out.versions.first().ifEmpty(null))
        }
    }

    /*
    ================================================================================================
                                                Binning
    ================================================================================================
    */

    ch_bowtie2_assembly_multiqc = Channel.empty()
    ch_busco_summary            = Channel.empty()
    ch_busco_multiqc            = Channel.empty()
    if (!params.skip_binning) {
        if (!params.single_end && !params.skip_spades) {
            METABAT2_BINNING_TRIMMOMATIC_SPADES (
            ch_spades_assemblies_trimmomatic,
            ch_short_reads_trimmomatic
        )
        ch_bowtie2_assembly_multiqc = ch_bowtie2_assembly_multiqc.mix(METABAT2_BINNING_TRIMMOMATIC_SPADES.out.bowtie2_assembly_multiqc)
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING_TRIMMOMATIC_SPADES.out.bowtie2_versions.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING_TRIMMOMATIC_SPADES.out.metabat2_versions.first().ifEmpty(null))
        }
        if (!params.skip_megahit) {
            METABAT2_BINNING_TRIMMOMATIC_MEGAHIT (
            ch_megahit_assemblies_trimmomatic,
            ch_short_reads_trimmomatic
        )
        ch_bowtie2_assembly_multiqc = ch_bowtie2_assembly_multiqc.mix(METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.bowtie2_assembly_multiqc)
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.bowtie2_versions.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.metabat2_versions.first().ifEmpty(null))
        }

        /*
         * BUSCO subworkflow: Quantitative measures for the assessment of genome assembly
         */
        if (!params.skip_busco) {
            if (!params.single_end && !params.skip_spades) {
                BUSCO_QC_TRIMMOMATIC_SPADES (
                    ch_busco_db_file,
                    ch_busco_download_folder,
                    METABAT2_BINNING_TRIMMOMATIC_SPADES.out.bins.transpose()
                )
                ch_busco_summary_trimmomatic_spades = BUSCO_QC_TRIMMOMATIC_SPADES.out.summary
                ch_busco_summary = ch_busco_summary.mix(ch_busco_summary_trimmomatic_spades)
                ch_busco_multiqc = ch_busco_multiqc.mix(BUSCO_QC_TRIMMOMATIC_SPADES.out.multiqc)
                ch_software_versions = ch_software_versions.mix(BUSCO_QC_TRIMMOMATIC_SPADES.out.versions.first().ifEmpty(null))

                // process information if BUSCO analysis failed for individual bins due to no matching genes
                BUSCO_QC_TRIMMOMATIC_SPADES.out
                    .failed_bin
                    .splitCsv(sep: '\t')
                    .map { bin, error -> if (!bin.contains(".unbinned.")) busco_failed_bins[bin] = error }
            }
            if (!params.skip_megahit) {
                BUSCO_QC_TRIMMOMATIC_MEGAHIT (
                    ch_busco_db_file,
                    ch_busco_download_folder,
                    METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.bins.transpose()
                )
                ch_busco_summary_trimmomatic_megahit = BUSCO_QC_TRIMMOMATIC_MEGAHIT.out.summary
                ch_busco_summary = ch_busco_summary.mix(ch_busco_summary_trimmomatic_megahit)
                ch_busco_multiqc = ch_busco_multiqc.mix(BUSCO_QC_TRIMMOMATIC_MEGAHIT.out.multiqc)
                ch_software_versions = ch_software_versions.mix(BUSCO_QC_TRIMMOMATIC_MEGAHIT.out.versions.first().ifEmpty(null))

                // process information if BUSCO analysis failed for individual bins due to no matching genes
                BUSCO_QC_TRIMMOMATIC_MEGAHIT.out
                    .failed_bin
                    .splitCsv(sep: '\t')
                    .map { bin, error -> if (!bin.contains(".unbinned.")) busco_failed_bins[bin] = error }
            }
        }

        if (!params.skip_quast){
            if (!params.single_end && !params.skip_spades) {
                QUAST_BINS_TRIMMOMATIC_SPADES ( METABAT2_BINNING_TRIMMOMATIC_SPADES.out.bins )
                ch_software_versions = ch_software_versions.mix(QUAST_BINS_TRIMMOMATIC_SPADES.out.versions.first().ifEmpty(null))

                QUAST_BINS_SUMMARY_TRIMMOMATIC_SPADES ( QUAST_BINS_TRIMMOMATIC_SPADES.out.quast_bin_summaries.collect() )
                ch_quast_bins_summary_trimmomatic_spades = QUAST_BINS_SUMMARY_TRIMMOMATIC_SPADES.out.summary
            }
            if (!params.skip_megahit) {
                QUAST_BINS_TRIMMOMATIC_MEGAHIT ( METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.bins )
                ch_software_versions = ch_software_versions.mix(QUAST_BINS_TRIMMOMATIC_MEGAHIT.out.versions.first().ifEmpty(null))

                QUAST_BINS_SUMMARY_TRIMMOMATIC_MEGAHIT ( QUAST_BINS_TRIMMOMATIC_MEGAHIT.out.quast_bin_summaries.collect() )
                ch_quast_bins_summary_trimmomatic_megahit = QUAST_BINS_SUMMARY_TRIMMOMATIC_MEGAHIT.out.summary
            }
        }

        /*
        -------------------------------------------------------------------------------------------------
                                Sourmash added for Taxonomic Classification
        -------------------------------------------------------------------------------------------------
        */

        /*
         * Sourmash: k-mer based taxonomic exploration and classification routines for genome and metagenome analysis
         */
        if (params.sourmash) {
            /* if (!params.single_end && !params.skip_spades) {
                SOURMASH_SIGNATURE_TRIMMOMATIC_SPADES ( METABAT2_BINNING_TRIMMOMATIC_SPADES.out.bins )

                SOURMASH_SUMMARIZE_TRIMMOMATIC_SPADES ( 
                    ch_sourmash_db,
                    SOURMASH_SIGNATURE_TRIMMOMATIC_SPADES.out.signatures
                )
                ch_software_versions = ch_software_versions.mix(SOURMASH_SIGNATURE_TRIMMOMATIC_SPADES.out.versions.first().ifEmpty(null))
            }
            if (!params.skip_megahit) {
                SOURMASH_SIGNATURE_TRIMMOMATIC_MEGAHIT ( METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.bins )

                SOURMASH_SUMMARIZE_TRIMMOMATIC_MEGAHIT ( 
                    ch_sourmash_db,
                    SOURMASH_SIGNATURE_TRIMMOMATIC_MEGAHIT.out.signatures
                )
                ch_software_versions = ch_software_versions.mix(SOURMASH_SIGNATURE_TRIMMOMATIC_MEGAHIT.out.versions.first().ifEmpty(null))
            } */

            if (!params.single_end && !params.skip_spades) {
                SOURMASH_TRIMMOMATIC_SPADES (ch_sourmash_db, METABAT2_BINNING_TRIMMOMATIC_SPADES.out.bins)
                ch_software_versions = ch_software_versions.mix(SOURMASH_TRIMMOMATIC_SPADES.out.versions.first().ifEmpty(null))
            }
            if (!params.skip_megahit) {
                SOURMASH_TRIMMOMATIC_MEGAHIT ( ch_sourmash_db, METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.bins )
                ch_software_versions = ch_software_versions.mix(SOURMASH_TRIMMOMATIC_MEGAHIT.out.versions.first().ifEmpty(null))
            }
        }

        /*
         * GTDB-tk: taxonomic classifications using GTDB reference
         */
        if (params.gtdbtk && gtdb) {
            ch_gtdbtk_summary_trimmomatic_spades = Channel.empty()
            ch_gtdbtk_summary_trimmomatic_megahit = Channel.empty()
            if (!params.single_end && !params.skip_spades) {
                GTDBTK_TRIMMOMATIC_SPADES ( 
                    METABAT2_BINNING_TRIMMOMATIC_SPADES.out.bins,
                    ch_busco_summary_trimmomatic_spades,
                    ch_gtdb
                )
                ch_software_versions = ch_software_versions.mix(GTDBTK_TRIMMOMATIC_SPADES.out.versions.first().ifEmpty(null))
                ch_gtdbtk_summary_trimmomatic_spades = GTDBTK_TRIMMOMATIC_SPADES.out.summary
                if (!params.skip_busco || !params.skip_quast || gtdb) {
                    BIN_SUMMARY_TRIMMOMATIC_SPADES (
                        METABAT2_BINNING_TRIMMOMATIC_SPADES.out.depths_summary,
                        ch_busco_summary_trimmomatic_spades.ifEmpty([]),
                        ch_quast_bins_summary_trimmomatic_spades.ifEmpty([]),
                        ch_gtdbtk_summary_trimmomatic_spades.ifEmpty([])
                    )
                    ch_bin_depths_summary = METABAT2_BINNING_TRIMMOMATIC_SPADES.out.depths_summary
                    ch_busco_summary = ch_busco_summary_trimmomatic_spades
                    ch_quast_bins_summary = ch_quast_bins_summary_trimmomatic_spades
                }
            }
            if (!params.skip_megahit) {
                GTDBTK_TRIMMOMATIC_MEGAHIT ( 
                    METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.bins,
                    ch_busco_summary_trimmomatic_megahit,
                    ch_gtdb
                )
                ch_software_versions = ch_software_versions.mix(GTDBTK_TRIMMOMATIC_MEGAHIT.out.versions.first().ifEmpty(null))
                ch_gtdbtk_summary_trimmomatic_megahit = GTDBTK_TRIMMOMATIC_MEGAHIT.out.summary
                if (!params.skip_busco || !params.skip_quast || gtdb) {
                    BIN_SUMMARY_TRIMMOMATIC_MEGAHIT (
                        METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.depths_summary,
                        ch_busco_summary_trimmomatic_megahit.ifEmpty([]),
                        ch_quast_bins_summary_trimmomatic_megahit.ifEmpty([]),
                        ch_gtdbtk_summary_trimmomatic_megahit.ifEmpty([])
                    )
                    ch_bin_depths_summary = METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.depths_summary
                    ch_busco_summary = ch_busco_summary_trimmomatic_megahit
                    ch_quast_bins_summary = ch_quast_bins_summary_trimmomatic_megahit
                }
            }
        }

        /*
         * Prokka: Genome annotation
         */
        if (!params.single_end && !params.skip_spades) {
            METABAT2_BINNING_TRIMMOMATIC_SPADES.out.bins
                .transpose()
                    .map { meta, bin ->
                        def meta_new = meta.clone()
                        meta_new.id  = bin.getBaseName()
                        [ meta_new, bin ]
                    }
                    .set { ch_bins_for_prokka_trimmomatic_spades }

            if (!params.skip_prokka){
                PROKKA_BINS_TRIMMOMATIC_SPADES (
                    ch_bins_for_prokka_trimmomatic_spades,
                    [],
                    []
                )
                PROKKA_ASSEMBLIES_TRIMMOMATIC_SPADES (
                    ch_spades_assemblies_trimmomatic,
                    [],
                    []
                )
                ch_software_versions = ch_software_versions.mix(PROKKA_BINS_TRIMMOMATIC_SPADES.out.versions.first().ifEmpty(null))
            }
        }
        if (!params.skip_megahit) {
            METABAT2_BINNING_TRIMMOMATIC_MEGAHIT.out.bins
                .transpose()
                    .map { meta, bin ->
                        def meta_new = meta.clone()
                        meta_new.id  = bin.getBaseName()
                        [ meta_new, bin ]
                    }
                    .set { ch_bins_for_prokka_trimmomatic_megahit }

            if (!params.skip_prokka){
                PROKKA_BINS_TRIMMOMATIC_MEGAHIT (
                    ch_bins_for_prokka_trimmomatic_megahit,
                    [],
                    []
                )
                PROKKA_ASSEMBLIES_TRIMMOMATIC_MEGAHIT (
                    ch_megahit_assemblies_trimmomatic,
                    [],
                    []
                )
                ch_software_versions = ch_software_versions.mix(PROKKA_BINS_TRIMMOMATIC_MEGAHIT.out.versions.first().ifEmpty(null))
            }
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    /* ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS_TRIMMOMATIC (
        ch_software_versions.map { it }.collect()
    ) */

    CUSTOM_DUMPSOFTWAREVERSIONS_TRIMMOMATIC (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMag.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS_TRIMMOMATIC.out.mqc_yml.collect())

    MULTIQC_TRIMMOMATIC (
        ch_multiqc_files.collect(),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIMMED_TRIMMOMATIC.out.zip.collect{it[1]}.ifEmpty([]),
        ch_bowtie2_removal_host_multiqc.collect{it[1]}.ifEmpty([]),
        ch_quast_multiqc.collect().ifEmpty([]),
        ch_bowtie2_assembly_multiqc.collect().ifEmpty([]),
        ch_busco_multiqc.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC_TRIMMOMATIC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC_TRIMMOMATIC.out.versions)

    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
   if (params.email || params.email_on_fail) {
       NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, busco_failed_bins)
  }
    NfcoreTemplate.summary(workflow, params, log, busco_failed_bins)
}

/*
========================================================================================
    THE END
========================================================================================
*/
