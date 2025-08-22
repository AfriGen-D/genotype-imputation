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
params.generate_chunk_plots = false
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
    ╔═══════════════════════════════════════════════════════════════════════════════╗
    ║                    ChiPImputation Genotype Imputation Pipeline                   ║
    ║                              Version ${params.version}                                       ║
    ╚═══════════════════════════════════════════════════════════════════════════════╝
    
    USAGE:
      nextflow run ${workflow.projectDir}/main_nfcore.nf --input samplesheet.csv --outdir results
    
    DESCRIPTION:
      A comprehensive Nextflow pipeline for genotype imputation with advanced QC,
      phasing, imputation, and reporting capabilities. Supports multiple reference
      panels and provides detailed quality metrics.
    
    ┌───────────────────────────────────────────────────────────────────────────────┐
    │ MANDATORY ARGUMENTS                                                          │
    └───────────────────────────────────────────────────────────────────────────────┘
      --input <file>        Path to input samplesheet CSV file
                           Format: dataset,vcf,population,study
      --outdir <path>       Directory where results will be saved
      
    ┌───────────────────────────────────────────────────────────────────────────────┐
    │ PROCESSING OPTIONS                                                           │
    └───────────────────────────────────────────────────────────────────────────────┘
      --chromosomes <str>   Chromosomes to process (default: ALL)
                           Options: ALL, 1-22, X, Y, or comma-separated list
      --chunk_size <int>    Size of genomic chunks in bp (default: 5000000)
      --buffer_size <int>   Buffer size for chunk boundaries (default: 500000)
      
    ┌───────────────────────────────────────────────────────────────────────────────┐
    │ QUALITY CONTROL                                                              │
    └───────────────────────────────────────────────────────────────────────────────┘
      --site_miss <float>   Max site missingness rate (default: 0.05)
      --hwe <float>         Hardy-Weinberg p-value threshold (default: 1e-5)
      --mac <int>           Minimum allele count (default: 1)
      --maf_thresh <float>  Minor allele frequency threshold (default: 0.1)
      --r2_threshold        R² threshold for quality (default: 0.3)
      
    ┌───────────────────────────────────────────────────────────────────────────────┐
    │ PHASING & IMPUTATION                                                         │
    └───────────────────────────────────────────────────────────────────────────────┘
      --phasing_method      Phasing algorithm: eagle, shapeit4 (default: eagle)
      --impute_method       Imputation: minimac4, impute5 (default: minimac4)
      --NE <int>            Effective population size (default: 20000)
      
    ┌───────────────────────────────────────────────────────────────────────────────┐
    │ EXECUTION PROFILES                                                           │
    └───────────────────────────────────────────────────────────────────────────────┘
      -profile docker       Use Docker containers
      -profile singularity  Use Singularity containers  
      -profile slurm        Submit jobs to SLURM cluster
      -profile test         Run with test dataset
      
    ┌───────────────────────────────────────────────────────────────────────────────┐
    │ EXAMPLES                                                                      │
    └───────────────────────────────────────────────────────────────────────────────┘
      # Basic run with singularity
      nextflow run main_nfcore.nf --input samples.csv --outdir results -profile singularity
      
      # SLURM cluster with specific chromosomes
      nextflow run main_nfcore.nf --input samples.csv --outdir results --chromosomes 1,2,3 -profile slurm,singularity
      
    ┌───────────────────────────────────────────────────────────────────────────────┐
    │ DOCUMENTATION & SUPPORT                                                      │
    └───────────────────────────────────────────────────────────────────────────────┘
      Documentation: https://github.com/h3abionet/chipimputation
      Issues:        https://github.com/h3abionet/chipimputation/issues
      Slack:         h3abionet.slack.com #imputation
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

// Determine executor configuration
def executorInfo = "local"
def queueInfo = "N/A"
def queueSizeInfo = "N/A"
if (workflow.profile.contains('slurm')) {
    executorInfo = "SLURM"
    queueInfo = workflow.configFiles.any { it.text.contains('process.queue') } ? "Main" : "default"
    queueSizeInfo = "100 jobs"
}

// ASCII art banner
def banner = """
╔═══════════════════════════════════════════════════════════════════════════════╗
║                                                                               ║
║    ░█████╗░██╗░░██╗██╗██████╗░██╗███╗░░░███╗██████╗░██╗░░░██╗████████╗        ║
║    ██╔══██╗██║░░██║██║██╔══██╗██║████╗░████║██╔══██╗██║░░░██║╚══██╔══╝        ║
║    ██║░░╚═╝███████║██║██████╔╝██║██╔████╔██║██████╔╝██║░░░██║░░░██║░░░        ║
║    ██║░░██╗██╔══██║██║██╔═══╝░██║██║╚██╔╝██║██╔═══╝░██║░░░██║░░░██║░░░        ║
║    ╚█████╔╝██║░░██║██║██║░░░░░██║██║░╚═╝░██║██║░░░░░╚██████╔╝░░░██║░░░        ║
║    ░╚════╝░╚═╝░░╚═╝╚═╝╚═╝░░░░░╚═╝╚═╝░░░░░╚═╝╚═╝░░░░░░╚═════╝░░░░╚═╝░░░        ║
║                                                                               ║
║         G E N O T Y P E   I M P U T A T I O N   P I P E L I N E               ║
║                          Version ${params.version}                            ║
╚═══════════════════════════════════════════════════════════════════════════════╝
""".stripIndent()

log.info banner

// Determine resource usage
def memoryGB = params.max_memory.toString().replaceAll(' GB', '')
def isHighMem = memoryGB.toInteger() >= 32
def performanceMode = isHighMem ? 'High Performance' : 'Standard'

// Calculate estimated variants
def estimatedVariants = params.chunk_size > 0 ? "~${(50000000 / params.chunk_size).round()} chunks per chromosome" : "N/A"

log.info """
┌───────────────────────────────────────────────────────────────────────────────┐
│ EXECUTION ENVIRONMENT                                                         │
├───────────────────────────────────────────────────────────────────────────────┤
│ Executor       : ${executorInfo.padRight(20)} Profile: ${workflow.profile ?: 'standard'}         │
│ Queue          : ${queueInfo.padRight(20)} Jobs:    ${queueSizeInfo}                   │
│ Work Directory : ${workflow.workDir}                                              
│ Mode           : ${performanceMode}                                               
└───────────────────────────────────────────────────────────────────────────────┘

┌───────────────────────────────────────────────────────────────────────────────┐
│ INPUT DATA                                                                    │
├───────────────────────────────────────────────────────────────────────────────┤
│ Samplesheet    : ${params.input}
│ Output Dir     : ${params.outdir}
│ Project        : ${params.project_name}
│ Genome Build   : ${params.genome_build ?: 'b38'}
└───────────────────────────────────────────────────────────────────────────────┘

┌───────────────────────────────────────────────────────────────────────────────┐
│ PROCESSING STRATEGY                                                           │
├───────────────────────────────────────────────────────────────────────────────┤
│ Chromosomes    : ${params.chromosomes == 'ALL' ? 'All autosomes + X' : params.chromosomes}
│ Chunking       : ${params.chunk_size} bp chunks with ${params.buffer_size} bp buffer
│ Estimated      : ${estimatedVariants}
└───────────────────────────────────────────────────────────────────────────────┘

┌───────────────────────────────────────────────────────────────────────────────┐
│ REFERENCE PANELS                                                              │
├───────────────────────────────────────────────────────────────────────────────┤
│ Panel Name     : ${params.ref_panels ? params.ref_panels[0][0] : 'None configured'}
│ Genetic Map    : ${params.eagle_genetic_map ? params.eagle_genetic_map.split('/')[-1] : 'Not specified'}
│ Reference      : ${params.reference_genome ? params.reference_genome.split('/')[-1] : 'Not specified'}
└───────────────────────────────────────────────────────────────────────────────┘

┌───────────────────────────────────────────────────────────────────────────────┐
│ QUALITY CONTROL THRESHOLDS                                                    │
├───────────────────────────────────────────────────────────────────────────────┤
│ Site Missingness : ≤ ${params.site_miss}        Hardy-Weinberg : p > ${params.hwe}
│ Min Allele Count : ≥ ${params.mac}             Min Alt Count  : ≥ ${params.min_ac}
│ MAF Threshold    : ≥ ${params.maf_thresh}        Max Mismatch   : ≤ ${params.max_mismatch_rate * 100}%
│ R² Threshold     : ≥ ${params.r2_threshold}        Info Cutoff    : ≥ ${params.impute_info_cutoff}
└───────────────────────────────────────────────────────────────────────────────┘

┌───────────────────────────────────────────────────────────────────────────────┐
│ ALGORITHMS                                                                    │
├───────────────────────────────────────────────────────────────────────────────┤
│ Phasing        : ${params.phasing_method.toUpperCase()} (${params.eagle_pbwt_iters} PBWT iterations)
│ Imputation     : ${params.impute_method.toUpperCase()} (Ne=${params.NE}, ${params.impute_iter} iterations, ${params.impute_burnin} burn-in)
│ Reporting      : ${params.aggregate_reports ? 'Hierarchical aggregation enabled' : 'Chunk-level only'}
│ Window Analysis: ${params.r2_window_size / 1000000} Mb windows for R² analysis
└───────────────────────────────────────────────────────────────────────────────┘

┌───────────────────────────────────────────────────────────────────────────────┐
│ COMPUTATIONAL RESOURCES                                                       │
├───────────────────────────────────────────────────────────────────────────────┤
│ Max CPUs       : ${params.max_cpus} cores
│ Max Memory     : ${params.max_memory}
│ Max Time       : ${params.max_time}
│ Error Strategy : ${params.max_retries} retries on failure
└───────────────────────────────────────────────────────────────────────────────┘

Starting pipeline execution...
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
        .filter { row ->
            // Skip rows where dataset starts with # (comment lines)
            !row.dataset.startsWith('#')
        }
        .map { row ->
            def meta = [:]
            meta.id = row.dataset
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
    // Call CHIPIMPUTATION directly to reduce nesting level
    CHIPIMPUTATION ( ch_input )
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