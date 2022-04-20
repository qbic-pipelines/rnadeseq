#!/usr/bin/env nextflow
/*
========================================================================================
    qbic-pipelines/rnadeseq
========================================================================================
    Github : https://github.com/qbic-pipelines/rnadeseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
print "lala"
print params.gtf
params.gtf = "hihihi"
print params.gtf
print params.organism
print WorkflowMain.getGenomeAttribute(params, 'gtf')
print WorkflowMain.getGenomeAttribute(params, 'organism')
print "schalalala"
if (!params.genome && (!params.gene_counts || !params.library || !params.gtf || !params.keytype)) { exit 1, 'Please provide either genome parameter or parameters for organism, library, gtf and keytype!' }
if (!params.organism) { params.organism = WorkflowMain.getGenomeAttribute(params, 'organism') }
if (!params.library) { params.library = WorkflowMain.getGenomeAttribute(params, 'library') }
if (!params.gtf) { params.gtf = WorkflowMain.getGenomeAttribute(params, 'gtf') }
if (!params.keytype) { params.keytype = WorkflowMain.getGenomeAttribute(params, 'keytype') }
print params.gtf
print params.organism
print "brö"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
WorkflowMain.initialise(workflow, params, log)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { RNADESEQ } from './workflows/rnadeseq'
//
// WORKFLOW: Run main qbic-pipelines/rnadeseq analysis pipeline
//
workflow QBIC_RNADESEQ {
    RNADESEQ ()
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    QBIC_RNADESEQ ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
