/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnadeseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.rawcounts, params.metadata, params.model,
    params.project_summary, params.versions,
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
// TODO: Channel names:
//ch_metadata_file_for_deseq2; ch_metadata_file_for_pathway = ch_metadata_file
//ch_model_for_deseq2_file; ch_model_for_report_file; ch_model_file_for_pathway = ch_model_file
//ch_genes_for_deseq2_file; ch_genes_for_report_file; ch_genes_for_pathway = ch_genes


// Check mandatory parameters
//TODO: Is Channel and channel the same?
if (params.rawcounts) { ch_counts_file = Channel.fromPath(params.rawcounts) } else { exit 1, 'Please provide raw counts file!' }
if (params.metadata) { ch_metadata_file = Channel.fromPath(params.metadata) } else { exit 1, 'Please provide metadata file!' }
if (params.model) { ch_model_file = Channel.fromPath(params.model) } else { exit 1, 'Please provide linear model file!' }
if (params.project_summary) { ch_proj_summary_file = Channel.fromPath(params.project_summary) } else { exit 1, 'Please provide project summary file!' }
if (params.versions) { ch_softwareversions_file = Channel.fromPath(params.versions) } else { exit 1, 'Please provide software versions file!' }


//TODO: is it necessary to add if (params.bla) etc. before each of these channels? The old pipeline did not check:
// Create channels for optional parameters
ch_contrast_matrix_for_deseq2 = Channel.fromPath(params.contrast_matrix)
ch_contrast_list_for_deseq2 = Channel.fromPath(params.contrast_list)
ch_contrast_pairs_for_deseq2 = Channel.fromPath(params.contrast_pairs)
ch_relevel_for_deseq2 = Channel.fromPath(params.relevel)
ch_quote_file = Channel.fromPath(params.quote)
ch_genes = Channel.fromPath(params.genelist)
ch_report_options_file = Channel.fromPath(params.report_options)
ch_kegg_blacklist_for_pathway = Channel.fromPath(params.kegg_blacklist)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
//TODO: Are any configs needed?
//ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
// TODO: Import here the four processes/modules
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RNADESEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
  //  INPUT_CHECK (
  //      ch_input
  //  )
  //  ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    // TODO: Here, I have to call the modules
    // This channel contains the versions of all tools of the current module

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}


/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
