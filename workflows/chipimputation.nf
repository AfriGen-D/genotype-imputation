/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESS            } from '../subworkflows/local/preprocess'
include { PHASE                 } from '../subworkflows/local/phase'
include { IMPUTE                } from '../subworkflows/local/impute'
include { REPORT                } from '../subworkflows/local/report'
include { REPORT_AGG   } from '../subworkflows/local/report_hierarchical'
// include { QC_PLOTS              } from '../subworkflows/local/qc_plots'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CHIPIMPUTATION {

    take:
    ch_input // channel: samplesheet read in from --input

    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // SUBWORKFLOW: Preprocessing and QC
    //
    PREPROCESS(
        ch_input
    )
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)
    
    //
    // SUBWORKFLOW: Phasing
    //
    PHASE(
        PREPROCESS.out.vcf
    )
    ch_versions = ch_versions.mix(PHASE.out.versions)
    
    //
    // SUBWORKFLOW: Imputation
    //
    IMPUTE(
        PHASE.out.phased
    )
    ch_versions = ch_versions.mix(IMPUTE.out.versions)
    
    //
    // SUBWORKFLOW: Reporting
    //
    // Combine imputed VCFs with their info files for reporting
    // Now passing full VCF files to enable frequency comparison plots
    ch_report_input = IMPUTE.out.imputed
        .join(IMPUTE.out.info, by: [0])  // Join by meta
        .map { meta, ref_name, vcf, vcf_index, info ->
            [meta, ref_name, vcf, vcf_index, info]  // Pass VCF files for frequency comparison
        }
    
    REPORT(
        ch_report_input,
        IMPUTE.out.ref_vcf
    )
    ch_versions = ch_versions.mix(REPORT.out.versions)
    
    //
    // SUBWORKFLOW: Hierarchical aggregation of reports (chunk → chromosome → genome)
    // This implements the three-tier aggregation system for post-imputation data
    //
    if (params.aggregate_reports != false) {
        // Collect just the well-imputed reports which have a consistent structure
        ch_chunk_reports = REPORT.out.well_imputed
            .map { meta, ref_name, well_imputed_file, summary_file ->
                def updated_meta = [:]
                updated_meta.id = meta.id
                updated_meta.sample = meta.sample ?: (meta.id.contains('_chr') ? meta.id.split('_chr')[0] : meta.id)
                updated_meta.dataset = meta.sample ?: (meta.id.contains('_chr') ? meta.id.split('_chr')[0] : meta.id)
                updated_meta.population = meta.population
                updated_meta.study = meta.study
                
                // Extract chromosome from chunk ID
                if (meta.id.contains('_chr')) {
                    def chr_part = meta.id.split('_chr')[1]
                    updated_meta.chromosome = chr_part.split('_')[0]
                } else {
                    updated_meta.chromosome = 'unknown'
                }
                
                // Return the summary file which contains the key metrics
                tuple(updated_meta, ref_name, summary_file)
            }
        
        // Collect performance plots which have consistent structure
        ch_chunk_plots = REPORT.out.performance_plots
            .map { meta, ref_name, plot_file ->
                def updated_meta = [:]
                updated_meta.id = meta.id
                updated_meta.sample = meta.sample ?: (meta.id.contains('_chr') ? meta.id.split('_chr')[0] : meta.id)
                updated_meta.dataset = meta.sample ?: (meta.id.contains('_chr') ? meta.id.split('_chr')[0] : meta.id)
                updated_meta.population = meta.population
                updated_meta.study = meta.study
                
                if (meta.id.contains('_chr')) {
                    def chr_part = meta.id.split('_chr')[1]
                    updated_meta.chromosome = chr_part.split('_')[0]
                } else {
                    updated_meta.chromosome = 'unknown'
                }
                
                tuple(updated_meta, ref_name, plot_file)
            }
        
        // Get info files for aggregation
        ch_chunk_info = REPORT.out.filtered_info
            .map { meta, ref_name, filtered_info, acc_info ->
                def updated_meta = [:]
                updated_meta.id = meta.id
                updated_meta.sample = meta.sample ?: (meta.id.contains('_chr') ? meta.id.split('_chr')[0] : meta.id)
                updated_meta.dataset = meta.sample ?: (meta.id.contains('_chr') ? meta.id.split('_chr')[0] : meta.id)
                updated_meta.population = meta.population
                updated_meta.study = meta.study
                
                if (meta.id.contains('_chr')) {
                    def chr_part = meta.id.split('_chr')[1]
                    updated_meta.chromosome = chr_part.split('_')[0]
                } else {
                    updated_meta.chromosome = 'unknown'
                }
                
                // Use the acc.info file which has the quality metrics
                tuple(updated_meta, ref_name, acc_info)
            }
        
        REPORT_AGG(
            ch_chunk_reports,
            ch_chunk_plots,
            ch_chunk_info
        )
        ch_versions = ch_versions.mix(REPORT_AGG.out.versions)
    }
    
    //
    // SUBWORKFLOW: QC Plots (disabled for now - needs implementation)
    //
    // if (params.qc_plots) {
    //     QC_PLOTS(
    //         IMPUTE.out.imputed
    //     )
    //     ch_versions = ch_versions.mix(QC_PLOTS.out.versions)
    //     ch_multiqc_files = ch_multiqc_files.mix(QC_PLOTS.out.plots)
    // }

    emit:
    imputed          = IMPUTE.out.imputed           // channel: [ val(meta), path(vcf) ]
    reports          = REPORT.out.reports           // channel: [ val(meta), path(reports) ]
    versions         = ch_versions                  // channel: [ path(versions.yml) ]
    multiqc          = ch_multiqc_files            // channel: [ path(files) ]
    // Hierarchical report outputs (when aggregate_reports is enabled)
    // chr_reports      = params.aggregate_reports ? REPORT_HIERARCHICAL.out.chr_reports : Channel.empty()
    // genome_reports   = params.aggregate_reports ? REPORT_HIERARCHICAL.out.genome_reports : Channel.empty()
    // final_reports    = params.aggregate_reports ? REPORT_HIERARCHICAL.out.final_reports : Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/