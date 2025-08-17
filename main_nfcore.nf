#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    h3abionet/chipimputation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/h3abionet/chipimputation
    Website: https://h3abionet.org
    Slack  : https://h3africa.slack.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER DEFAULTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define default parameters
params.help = false
params.input = null
params.outdir = './results'
params.chromosomes = 'ALL'
params.chunk_size = 5000000
params.minRatio = 0.01
params.qc_plots = true
params.max_cpus = 16
params.max_memory = '128.GB'
params.version = '1.0.0'
params.project_name = 'h3achipimputation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Print help message if needed
if (params.help) {
    log.info """
    =========================================
    h3abionet/chipimputation v${params.version}
    =========================================
    
    Usage:
    nextflow run ${workflow.projectDir}/main_nfcore.nf --input samplesheet.csv --outdir results -profile docker
    
    Mandatory arguments:
      --input               Path to input samplesheet CSV file
      --outdir              The output directory where results will be saved
      
    Optional arguments:
      --chromosomes         Chromosomes to process (default: ALL)
      --chunk_size          Chunk size for processing (default: 5000000)
      --minRatio           Minimum ratio for imputation (default: 0.01)
      --qc_plots           Generate QC plots (default: true)
      
    Profiles:
      -profile docker       Use Docker containers
      -profile singularity  Use Singularity containers
      -profile conda        Use Conda environments
      -profile test         Run with test dataset
    """.stripIndent()
    System.exit(0)
}

// Check required parameters
if (!params.input) {
    log.error "Please provide an input samplesheet with --input"
    System.exit(1)
}

if (!params.outdir) {
    log.error "Please provide an output directory with --outdir"
    System.exit(1)
}

log.info """
=========================================
h3abionet/chipimputation v${params.version}
=========================================
Input Configuration:
  Input file     : ${params.input}
  Output dir     : ${params.outdir}
  Project name   : ${params.project_name}

Processing Parameters:
  Chromosomes    : ${params.chromosomes}
  Chunk size     : ${params.chunk_size} bp
  Buffer size    : ${params.buffer_size} bp

Reference Panels:
  Panel name     : ${params.ref_panels ? params.ref_panels[0][0] : 'None'}
  Genetic map    : ${params.eagle_genetic_map}
  Reference genome: ${params.reference_genome}

Quality Control:
  Site miss      : ${params.site_miss}
  HWE threshold  : ${params.hwe}
  Min allele count: ${params.mac}
  Min alt count  : ${params.min_ac}
  MAF threshold  : ${params.maf_thresh}
  Max mismatch   : ${params.max_mismatch_rate * 100}%

Overlap Checking:
  Min ratio      : ${params.minRatio}
  Min overlap    : 50 variants (hardcoded)

Phasing:
  Method         : ${params.phasing_method}
  PBWT iterations: ${params.eagle_pbwt_iters}

Imputation:
  Method         : ${params.impute_method}
  NE (pop size)  : ${params.NE}
  Iterations     : ${params.impute_iter}
  Burn-in        : ${params.impute_burnin}
  Info cutoff    : ${params.impute_info_cutoff}
  R2 threshold   : ${params.r2_threshold}

Resources:
  Max CPUs       : ${params.max_cpus}
  Max memory     : ${params.max_memory}
  Max time       : ${params.max_time}
=========================================
"""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE INPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Create input channel from samplesheet CSV
if (params.input) {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.population = row.population ?: 'ALL'
            meta.study = row.study ?: params.project_name
            [ meta, file(row.vcf) ]
        }
        .set { ch_input }
} else {
    error "Please provide an input samplesheet with --input"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHIPIMPUTATION } from './workflows/chipimputation'

//
// WORKFLOW: Run main h3abionet/chipimputation analysis pipeline
//
workflow H3ABIONET_CHIPIMPUTATION {
    
    take:
    ch_input  // channel: [ val(meta), path(vcf) ]
    
    main:
    CHIPIMPUTATION ( ch_input )
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
    H3ABIONET_CHIPIMPUTATION ( ch_input )
}

//
// Print summary of skipped chunks on completion
//
workflow.onComplete {
    log.info """
    ========================================
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    ========================================
    
    Check the reports directory for detailed information about:
    - Chunks skipped due to insufficient overlap
    - Chunks skipped due to high allele mismatch
    - Overall imputation quality metrics
    
    Reports location: ${params.outdir}/reports/
    ========================================
    """.stripIndent()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/