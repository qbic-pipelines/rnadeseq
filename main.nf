#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    qbic-pipelines/rnadeseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/qbic-pipelines/rnadeseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_rnadeseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//Either genome needs to be set or the parameters gtf (for rsem/salmon) and organism, species_library and keytype
//(for pathway analysis) have to be provided separately
//-->gtf is necessary for rsem and salmon
if (params.input_type in ["rsem", "salmon"]) {
    if (!params.genome && !params.gtf) { exit 1, 'Please provide either genome or gtf file!' }
    else if (!params.gtf) {
        params.gtf = getGenomeAttribute('gtf')
        if (!params.gtf) {
            exit 1, 'It seems that for your genome, no gtf file is defined. Please provide a gtf with the --gtf parameter or open a github issue: https://github.com/qbic-pipelines/rnadeseq/issues'
        }
    }
}

//-->organism, species_library and keytype are necessary for pathway analysis
if (params.run_pathway_analysis) {
    if (!params.genome && !params.organism) { exit 1, 'Please provide either genome or organism!' }
    else if (!params.organism) {
        params.organism = getGenomeAttribute('organism')
        if (!params.organism) {
            exit 1, 'It seems that for your genome, no organism is defined. Please provide the organism with the --organism parameter or open a github issue: https://github.com/qbic-pipelines/rnadeseq/issues'
        }
    }
    if (!params.genome && !params.species_library) { exit 1, 'Please provide either genome or species_library!' }
    else if (!params.species_library) {
        params.species_library = getGenomeAttribute('species_library')
        if (!params.species_library) {
            exit 1, 'It seems that for your genome, no species_library is defined. Please provide the library with the --species_library parameter or open a github issue: https://github.com/qbic-pipelines/rnadeseq/issues'
        }
    }
    if (!params.genome && !params.keytype) { exit 1, 'Please provide either genome or keytype!' }
    else if (!params.keytype) {
        params.keytype = getGenomeAttribute('keytype')
        if (!params.keytype) {
            exit 1, 'It seems that for your genome, no keytype is defined. Please provide the keytype with the --keytype parameter or open a github issue: https://github.com/qbic-pipelines/rnadeseq/issues'
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNADESEQ } from './workflows/rnadeseq'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow QBIC_RNADESEQ {

    main:

    //
    // WORKFLOW: Run pipeline
    //
    RNADESEQ ()

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // WORKFLOW: Run main workflow
    //
    QBIC_RNADESEQ ()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
