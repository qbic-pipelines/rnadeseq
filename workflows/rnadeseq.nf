/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.model
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (!(params.input_type in ["featurecounts", "salmon", "rsem", "smrnaseq"])) { exit 1, 'Wrong input type ' + params.input_type + ', must be one of "featurecounts", "salmon", "rsem", "smrnaseq"!' }

ch_counts_path = Channel.fromPath(params.gene_counts)
ch_metadata_file = Channel.fromPath(params.input)
ch_model_file = Channel.fromPath(params.model)

// Create channel for genome parameter gtf (the other genome params are not files)
if (params.input_type in ["rsem", "salmon"]) { ch_gtf = Channel.fromPath(params.gtf) } else { ch_gtf = Channel.fromPath("FALSE") }

// Create channels for optional parameters
ch_contrast_matrix = Channel.fromPath(params.contrast_matrix)
ch_contrast_list = Channel.fromPath(params.contrast_list)
ch_contrast_pairs = Channel.fromPath(params.contrast_pairs)
ch_relevel = Channel.fromPath(params.relevel)
ch_genes = Channel.fromPath(params.genelist)
ch_multiqc_file = Channel.fromPath(params.multiqc)
ch_custom_gmt = Channel.fromPath(params.custom_gmt)
ch_custom_background = Channel.fromPath(params.custom_background)
ch_proj_summary_file = Channel.fromPath(params.project_summary)
ch_softwareversions_file = Channel.fromPath(params.software_versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REPORT } from '../modules/local/report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNADESEQ {

    main:

    REPORT (
        ch_counts_path,
        ch_metadata_file,
        ch_model_file,
        ch_gtf,

        ch_contrast_matrix,
        ch_contrast_list,
        ch_contrast_pairs,
        ch_genes,
        ch_relevel,

        ch_proj_summary_file,
        ch_softwareversions_file,
        ch_multiqc_file,
        ch_custom_gmt,
        ch_custom_background
    )

    //TODO: Enable this:
    // This channel contains the versions of all tools of the current module
//    CUSTOM_DUMPSOFTWAREVERSIONS (
//        ch_softwareversions_file.unique().collectFile(name: 'collated_versions.yml')
//    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
