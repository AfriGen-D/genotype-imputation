/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORT_BY_CHR
    Chromosome-level and genome-wide reporting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FILTER_INFO_BY_TARGET     } from '../../modules/local/report/filter_info_by_target'
include { REPORT_WELL_IMPUTED       } from '../../modules/local/report/report_well_imputed'
include { PLOT_PERFORMANCE          } from '../../modules/local/report/plot_performance'
include { REPORT_ACCURACY           } from '../../modules/local/report/report_accuracy'
include { PLOT_ACCURACY             } from '../../modules/local/report/plot_accuracy'
include { PLOT_R2_MAF               } from '../../modules/local/report/plot_r2_maf'
include { PLOT_R2_SNPPOS            } from '../../modules/local/report/plot_r2_snppos'
include { PLOT_R2_SNPCOUNT          } from '../../modules/local/report/plot_r2_snpcount'
include { PLOT_HIST_R2_SNPCOUNT     } from '../../modules/local/report/plot_hist_r2_snpcount'
include { PLOT_MAF_R2               } from '../../modules/local/report/plot_maf_r2'
include { AVERAGE_R2                } from '../../modules/local/report/average_r2'
include { PLOT_DOSAGE_DISTRIBUTION  } from '../../modules/local/report/plot_dosage_distribution'
include { PLOT_CALIBRATION          } from '../../modules/local/report/plot_calibration'
include { PLOT_CONCORDANCE_MAF      } from '../../modules/local/report/plot_concordance_maf'
include { PLOT_CROSS_VALIDATION     } from '../../modules/local/report/plot_cross_validation'
include { PLOT_HETEROZYGOSITY       } from '../../modules/local/report/plot_heterozygosity'
include { PLOT_HWE_DEVIATION        } from '../../modules/local/report/plot_hwe_deviation'
include { GENERATE_SUMMARY_REPORT   } from '../../modules/local/report/generate_summary'
include { COLLECT_WARNINGS          } from '../../modules/local/report/collect_warnings'

workflow REPORT_BY_CHR {
    take:
    ch_imputed_chr // channel: [ val(meta), val(ref_name), path(info), val(chr) ]
    ch_imputed_vcf // channel: [ val(meta), val(ref_name), path(vcf), val(chr) ]
    ch_true_vcf    // channel: [ val(meta), val(ref_name), path(true_vcf), val(chr) ] - optional

    main:
    ch_versions = Channel.empty()
    
    //
    // Split channels into chromosome-specific and genome-wide
    //
    ch_chr_data = ch_imputed_chr
        .filter { meta, ref_name, info, chr -> chr != null && chr != '' }
    
    ch_genome_data = ch_imputed_chr
        .filter { meta, ref_name, info, chr -> chr == null || chr == '' }
    
    ch_vcf_chr = ch_imputed_vcf
        .filter { meta, ref_name, vcf, chr -> chr != null && chr != '' }
    
    ch_vcf_genome = ch_imputed_vcf
        .filter { meta, ref_name, vcf, chr -> chr == null || chr == '' }
    
    //
    // CHROMOSOME-LEVEL ANALYSIS
    //
    
    // Basic R² metrics per chromosome
    PLOT_R2_SNPPOS ( ch_chr_data )
    ch_versions = ch_versions.mix(PLOT_R2_SNPPOS.out.versions)
    
    PLOT_R2_SNPCOUNT ( ch_chr_data )
    ch_versions = ch_versions.mix(PLOT_R2_SNPCOUNT.out.versions)
    
    PLOT_HIST_R2_SNPCOUNT ( ch_chr_data )
    ch_versions = ch_versions.mix(PLOT_HIST_R2_SNPCOUNT.out.versions)
    
    PLOT_MAF_R2 ( ch_chr_data )
    ch_versions = ch_versions.mix(PLOT_MAF_R2.out.versions)
    
    AVERAGE_R2 ( ch_chr_data )
    ch_versions = ch_versions.mix(AVERAGE_R2.out.versions)
    
    // Advanced QC plots per chromosome (if VCF available)
    if (ch_vcf_chr) {
        PLOT_DOSAGE_DISTRIBUTION ( ch_vcf_chr )
        ch_versions = ch_versions.mix(PLOT_DOSAGE_DISTRIBUTION.out.versions)
        
        PLOT_HETEROZYGOSITY ( ch_vcf_chr )
        ch_versions = ch_versions.mix(PLOT_HETEROZYGOSITY.out.versions)
        
        PLOT_HWE_DEVIATION ( ch_vcf_chr )
        ch_versions = ch_versions.mix(PLOT_HWE_DEVIATION.out.versions)
    }
    
    // Calibration per chromosome
    PLOT_CALIBRATION ( ch_chr_data )
    ch_versions = ch_versions.mix(PLOT_CALIBRATION.out.versions)
    
    // Concordance analysis per chromosome (if validation data available)
    ch_concordance = ch_vcf_chr
        .join(ch_true_vcf, by: [0, 1, 3], remainder: true)
        .map { meta, ref_name, chr, imputed_vcf, true_vcf ->
            [ meta, ref_name, imputed_vcf, true_vcf ?: file('NO_FILE'), chr ]
        }
    
    PLOT_CONCORDANCE_MAF ( ch_concordance )
    ch_versions = ch_versions.mix(PLOT_CONCORDANCE_MAF.out.versions)
    
    //
    // GENOME-WIDE ANALYSIS (aggregate across all chromosomes)
    //
    
    // Collect all chromosome data for genome-wide analysis
    ch_all_chr_info = ch_chr_data
        .map { meta, ref_name, info, chr -> [ meta, ref_name, info ] }
        .groupTuple(by: [0, 1])
        .map { meta, ref_name, info_files ->
            [ meta, ref_name, info_files, null ]  // null for genome-wide
        }
    
    // Genome-wide R² metrics
    PLOT_R2_SNPPOS.genome ( ch_all_chr_info )
    PLOT_R2_SNPCOUNT.genome ( ch_all_chr_info )
    PLOT_HIST_R2_SNPCOUNT.genome ( ch_all_chr_info )
    PLOT_MAF_R2.genome ( ch_all_chr_info )
    AVERAGE_R2.genome ( ch_all_chr_info )
    
    // Genome-wide advanced QC
    ch_all_chr_vcf = ch_vcf_chr
        .map { meta, ref_name, vcf, chr -> [ meta, ref_name, vcf ] }
        .groupTuple(by: [0, 1])
        .map { meta, ref_name, vcf_files ->
            // Merge VCFs or select representative
            [ meta, ref_name, vcf_files[0], null ]  // Using first VCF as representative
        }
    
    if (ch_all_chr_vcf) {
        PLOT_DOSAGE_DISTRIBUTION.genome ( ch_all_chr_vcf )
        PLOT_HETEROZYGOSITY.genome ( ch_all_chr_vcf )
        PLOT_HWE_DEVIATION.genome ( ch_all_chr_vcf )
    }
    
    PLOT_CALIBRATION.genome ( ch_all_chr_info )
    
    // Cross-validation plot (genome-wide only)
    ch_cv_results = Channel.value(file('NO_FILE'))  // Placeholder for CV results
    ch_cv_input = ch_all_chr_info
        .map { meta, ref_name, info, chr ->
            [ meta, ref_name, file('NO_FILE'), null ]
        }
    
    PLOT_CROSS_VALIDATION ( ch_cv_input )
    ch_versions = ch_versions.mix(PLOT_CROSS_VALIDATION.out.versions)
    
    //
    // GENERATE COMPREHENSIVE SUMMARY REPORT
    //
    
    // Collect all metrics and plots
    ch_all_metrics = AVERAGE_R2.out.average
        .mix(AVERAGE_R2.out.summary)
        .mix(AVERAGE_R2.genome.out.average)
        .mix(AVERAGE_R2.genome.out.summary)
        .groupTuple(by: [0, 1])
    
    ch_all_plots = PLOT_R2_SNPPOS.out.plot
        .mix(PLOT_R2_SNPCOUNT.out.plot)
        .mix(PLOT_HIST_R2_SNPCOUNT.out.plot)
        .mix(PLOT_MAF_R2.out.plot)
        .mix(PLOT_DOSAGE_DISTRIBUTION.out.plot)
        .mix(PLOT_CALIBRATION.out.plot)
        .mix(PLOT_CONCORDANCE_MAF.out.plot)
        .mix(PLOT_HETEROZYGOSITY.out.plot)
        .mix(PLOT_HWE_DEVIATION.out.plot)
        .mix(PLOT_CROSS_VALIDATION.out.plot)
        .mix(PLOT_R2_SNPPOS.genome.out.plot)
        .mix(PLOT_R2_SNPCOUNT.genome.out.plot)
        .mix(PLOT_HIST_R2_SNPCOUNT.genome.out.plot)
        .mix(PLOT_MAF_R2.genome.out.plot)
        .groupTuple(by: [0, 1])
    
    ch_summary_input = ch_all_metrics
        .join(ch_all_plots, by: [0, 1], remainder: true)
        .map { meta, ref_name, metrics, plots ->
            [ meta, ref_name, metrics ?: [file('NO_FILES')], plots ?: [file('NO_FILES')] ]
        }
    
    GENERATE_SUMMARY_REPORT ( ch_summary_input )
    ch_versions = ch_versions.mix(GENERATE_SUMMARY_REPORT.out.versions)
    
    emit:
    chr_reports     = AVERAGE_R2.out.average
    genome_reports  = AVERAGE_R2.genome.out.average
    chr_plots       = PLOT_R2_SNPPOS.out.plot
    genome_plots    = PLOT_R2_SNPPOS.genome.out.plot
    summary_report  = GENERATE_SUMMARY_REPORT.out.report
    versions        = ch_versions
}