#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
=                          Comprehensive QC Plotting Module                            =
========================================================================================
 Additional quality control plots for imputation assessment
----------------------------------------------------------------------------------------
*/

// Dosage Distribution Plot
process plot_dosage_distribution {
    tag "dosage_dist_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/dosage", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_file}")
    
    script:
        output_file = "${dataset}_${ref_panel}_dosage_distribution.pdf"
        template "dosage_distribution_plot.py"
}

// INFO Score Distribution Plot
process plot_info_distribution {
    tag "info_dist_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/info_score", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_file}")
    
    script:
        output_file = "${dataset}_${ref_panel}_info_distribution.pdf"
        template "info_score_distribution.py"
}

// Concordance vs MAF Plot
process plot_concordance_maf {
    tag "concordance_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/concordance", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file), file(masked_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_file}")
    
    script:
        output_file = "${dataset}_${ref_panel}_concordance_maf.pdf"
        masked_file = masked_file.name != 'NO_FILE' ? masked_file : ""
        template "concordance_maf_plot.py"
}

// Calibration Plot
process plot_calibration {
    tag "calibration_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/calibration", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_file}")
    
    script:
        output_file = "${dataset}_${ref_panel}_calibration.pdf"
        template "calibration_plot.py"
}

// Cross-validation Plot
process plot_cross_validation {
    tag "cross_val_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/cross_validation", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file), file(cv_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_file}")
    
    script:
        output_file = "${dataset}_${ref_panel}_cross_validation.pdf"
        cv_file = cv_file.name != 'NO_FILE' ? cv_file : ""
        template "cross_validation_plot.py"
}

// Heterozygosity Plot
process plot_heterozygosity {
    tag "heterozygosity_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/heterozygosity", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file), file(vcf_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_file}")
    
    script:
        output_file = "${dataset}_${ref_panel}_heterozygosity.pdf"
        vcf_file = vcf_file.name != 'NO_FILE' ? vcf_file : ""
        template "heterozygosity_plot.py"
}

// Hardy-Weinberg Equilibrium Plot
process plot_hwe_deviation {
    tag "hwe_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/hwe", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file), file(vcf_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_file}")
    
    script:
        output_file = "${dataset}_${ref_panel}_hwe_deviation.pdf"
        vcf_file = vcf_file.name != 'NO_FILE' ? vcf_file : ""
        template "hwe_deviation_plot.py"
}

// Ancestry-aware Assessment
process plot_ancestry_assessment {
    tag "ancestry_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/ancestry", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file), file(ancestry_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_prefix}_quality_assessment.pdf"), file("${output_prefix}_metrics.tsv")
    
    script:
        output_prefix = "${dataset}_${ref_panel}_ancestry"
        ancestry_file = ancestry_file.name != 'NO_FILE' ? ancestry_file : ""
        template "ancestry_imputation_assessment.py"
}

// Comprehensive QC Report (combines all plots)
process generate_comprehensive_qc {
    tag "comprehensive_qc_${dataset}_${ref_panel}"
    publishDir "${params.outDir}/qc_plots/comprehensive", overwrite: true, mode:'copy'
    label "python_plotting"
    label "bigmem"
    
    input:
        tuple val(dataset), val(ref_panel), file(info_file), file(masked_file), file(cv_file)
    
    output:
        tuple val(dataset), val(ref_panel), file("${output_prefix}_comprehensive_qc.pdf"), file("${output_prefix}_summary_stats.tsv")
    
    script:
        output_prefix = "${dataset}_${ref_panel}"
        masked_file = masked_file.name != 'NO_FILE' ? masked_file : ""
        cv_file = cv_file.name != 'NO_FILE' ? cv_file : ""
        template "comprehensive_imputation_qc.py"
}

// Workflow to run all QC plots
workflow run_qc_plots {
    take:
        imputation_data  // tuple of (dataset, ref_panel, info_file)
        masked_data      // Optional: masked variant data
        cv_data          // Optional: cross-validation data
        ancestry_data    // Optional: ancestry information
    
    main:
        // Basic QC plots
        plot_dosage_distribution(imputation_data)
        plot_info_distribution(imputation_data)
        plot_calibration(imputation_data)
        
        // Plots requiring additional data
        if (masked_data) {
            concordance_input = imputation_data.join(masked_data, by: [0, 1])
            plot_concordance_maf(concordance_input)
        }
        
        if (cv_data) {
            cv_input = imputation_data.join(cv_data, by: [0, 1])
            plot_cross_validation(cv_input)
        }
        
        if (ancestry_data) {
            ancestry_input = imputation_data.join(ancestry_data, by: [0, 1])
            plot_ancestry_assessment(ancestry_input)
        }
        
        // Heterozygosity and HWE (can use VCF if available)
        vcf_data = Channel.value(file('NO_FILE'))
        het_input = imputation_data.combine(vcf_data)
        plot_heterozygosity(het_input)
        plot_hwe_deviation(het_input)
        
        // Comprehensive report
        comprehensive_input = imputation_data
            .combine(masked_data ?: Channel.value(file('NO_FILE')))
            .combine(cv_data ?: Channel.value(file('NO_FILE')))
        generate_comprehensive_qc(comprehensive_input)
    
    emit:
        dosage_plots = plot_dosage_distribution.out
        info_plots = plot_info_distribution.out
        calibration_plots = plot_calibration.out
        comprehensive_report = generate_comprehensive_qc.out
}