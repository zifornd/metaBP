


include { SOURMASH_SIGNATURE             } from '../../modules/local/sourmash_signature'          //addParams( options: modules['sourmash_signature_cutadapt_spades']                    )
include { SOURMASH_SUMMARIZE             } from '../../modules/local/sourmash_summarize'          //addParams( options: modules['sourmash_summarize_cutadapt_spades']                    )

/* include { SOURMASH_SIGNATURE as SOURMASH_SIGNATURE_CUTADAPT_SPADES            } from '../modules/local/sourmash_signature'       addParams( options: modules['sourmash_signature_cutadapt_spades']                    )
include { SOURMASH_SIGNATURE as SOURMASH_SIGNATURE_CUTADAPT_MEGAHIT           } from '../modules/local/sourmash_signature'          addParams( options: modules['sourmash_signature_cutadapt_megahit']                   )
include { SOURMASH_SIGNATURE as SOURMASH_SIGNATURE_TRIMMOMATIC_SPADES         } from '../modules/local/sourmash_signature'          addParams( options: modules['sourmash_signature_trimmomatic_spades']                 )
include { SOURMASH_SIGNATURE as SOURMASH_SIGNATURE_TRIMMOMATIC_MEGAHIT        } from '../modules/local/sourmash_signature'          addParams( options: modules['sourmash_signature_trimmomatic_megahit']                )
include { SOURMASH_SUMMARIZE as SOURMASH_SUMMARIZE_CUTADAPT_SPADES            } from '../modules/local/sourmash_summarize'          addParams( options: modules['sourmash_summarize_cutadapt_spades']                    )
include { SOURMASH_SUMMARIZE as SOURMASH_SUMMARIZE_CUTADAPT_MEGAHIT           } from '../modules/local/sourmash_summarize'          addParams( options: modules['sourmash_summarize_cutadapt_megahit']                   )
include { SOURMASH_SUMMARIZE as SOURMASH_SUMMARIZE_TRIMMOMATIC_SPADES         } from '../modules/local/sourmash_summarize'          addParams( options: modules['sourmash_summarize_trimmomatic_spades']                 )
include { SOURMASH_SUMMARIZE as SOURMASH_SUMMARIZE_TRIMMOMATIC_MEGAHIT        } from '../modules/local/sourmash_summarize'          addParams( options: modules['sourmash_summarize_trimmomatic_megahit']                )
 */
workflow SOURMASH {
    take:
    ch_sourmash_db
    bins

    main:

    SOURMASH_SIGNATURE ( bins )
                
    SOURMASH_SUMMARIZE ( 
        ch_sourmash_db,
        SOURMASH_SIGNATURE.out.signatures
    )

    emit:
    versions = SOURMASH_SIGNATURE.out.versions
    sourmash_signature = SOURMASH_SIGNATURE.out.signatures
    sourmash_classification = SOURMASH_SUMMARIZE.out.report
}