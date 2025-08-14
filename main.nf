#!/usr/bin/env nextflow
/*
========================================================================================
    h3abionet/chipimputation
========================================================================================
    Github : https://github.com/h3abionet/chipimputation
    Website: https://h3africa.org
    Slack  : https://h3abionetwgs.slack.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

// Validate essential parameters
if (!params.input) {
    log.error "Please provide an input samplesheet with --input"
    System.exit(1)
}

if (!params.reference_panels || params.reference_panels.size() == 0) {
    log.error "Please provide reference panels with --reference_panels"
    System.exit(1)
}

// Print workflow header
log.info """
====================================
 GENOTYPE IMPUTATION PIPELINE v2.0
====================================
Input            : ${params.input}
Output directory : ${params.outdir}
Genome build     : ${params.genome_build}
Phasing tool     : ${params.phasing_tool}
Imputation tool  : ${params.imputation_tool}
Work directory   : ${workflow.workDir}
Profile          : ${workflow.profile}
====================================
""".stripIndent()

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { GENOTYPE_IMPUTATION } from './workflows/genotype_imputation'

//
// WORKFLOW: Run main h3abionet/chipimputation analysis pipeline
//
workflow H3ABIONET_CHIPIMPUTATION {
    GENOTYPE_IMPUTATION()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    H3ABIONET_CHIPIMPUTATION()
}

/*
========================================================================================
    THE END
========================================================================================
*/