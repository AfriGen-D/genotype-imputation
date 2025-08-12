#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
=                     Chromosome-level Reporting Module                                =
========================================================================================
 Chromosome-specific imputation quality reports and plots
----------------------------------------------------------------------------------------
*/

// Filter info files by chromosome
process filter_info_by_target_chr {
    tag "filter_${dataset_name}_${tagName}_${ref_panels.join('-')}_chr${chr}"
    label "bigmem"
    publishDir "${params.outDir}/reports_chr/chr${chr}/${ref_panels}", overwrite: true, mode:'copy'

    input:
        tuple val(dataset_name), val(ref_panels), val(chr), val(ref_infos)
    output:
        tuple val(dataset_name), val(ref_panels), val(chr), file("${well_out}_${info_cutoff}.tsv"), file("${acc_out}_${info_cutoff}.tsv")
    script:
        tagName = "${ref_panels.join('-')}"
        comb_info = "${dataset_name}_${tagName}_chr${chr}.imputed_info"
        well_out = "${comb_info}_well_imputed"
        acc_out = "${comb_info}_accuracy"
        infos = ref_infos  // Already joined
        datasets = ref_panels
        info_cutoff = params.impute_info_cutoff
        """
        python3 ${projectDir}/templates/improved/filter_info_minimac.py \\
            --infoFiles ${infos} \\
            --datasets ${datasets} \\
            --out_prefix ${comb_info} \\
            --infoCutoff ${info_cutoff}
        """
}

// Duplicate for dataset grouping (different tagging)
process filter_info_by_target_chr2 {
    tag "filter_${dataset_name}_${ref_panels.join('-')}_chr${chr}"
    label "bigmem"
    publishDir "${params.outDir}/reports_chr/chr${chr}/${dataset_name}", overwrite: true, mode:'copy'

    input:
        tuple val(dataset_name), val(ref_panels), val(chr), val(ref_infos)
    output:
        tuple val(dataset_name), val(ref_panels), val(chr), file("${well_out}_${info_cutoff}.tsv"), file("${acc_out}_${info_cutoff}.tsv")
    script:
        tagName = "${ref_panels.join('-')}"
        comb_info = "${dataset_name}_${tagName}_chr${chr}.imputed_info"
        well_out = "${comb_info}_well_imputed"
        acc_out = "${comb_info}_accuracy"
        infos = ref_infos  // Already joined
        datasets = ref_panels
        info_cutoff = params.impute_info_cutoff
        """
        python3 ${projectDir}/templates/improved/filter_info_minimac.py \\
            --infoFiles ${infos} \\
            --datasets ${datasets} \\
            --out_prefix ${comb_info} \\
            --infoCutoff ${info_cutoff}
        """
}

// Report well imputed by target - chromosome level
process report_well_imputed_by_target_chr {
    tag "report_wellImputed_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/reports_chr/chr${chr}/${ref_panels}", overwrite: true, mode:'copy'
    label "medium"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), file(inWell_imputed)
    
    output:
        tuple val(target_name), val(ref_panels), val(chr), file("${out_prefix}.tsv"), file("${out_prefix}_summary.tsv")
    
    script:
        out_prefix = "${inWell_imputed.baseName}.chr${chr}.imputed_info_performance_by_maf_report"
        datasets = ref_panels.split(',').join(',')
        """
        python3 ${projectDir}/templates/improved/report_well_imputed.py \\
            --info-files ${inWell_imputed} \\
            --datasets ${datasets} \\
            --output-prefix ${out_prefix} \\
            --rsq-threshold ${params.impute_info_cutoff} \\
            --maf-threshold ${params.maf_thresh}
        """
}

// Duplicate for dataset grouping
process report_well_imputed_by_target_chr2 {
    tag "report_wellImputed_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/reports_chr/chr${chr}/${target_name}", overwrite: true, mode:'copy'
    label "medium"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), file(inWell_imputed)
    
    output:
        tuple val(target_name), val(ref_panels), val(chr), file("${out_prefix}.tsv"), file("${out_prefix}_summary.tsv")
    
    script:
        out_prefix = "${inWell_imputed.baseName}.chr${chr}.imputed_info_performance_by_maf_report"
        datasets = ref_panels.split(',').join(',')
        """
        python3 ${projectDir}/templates/improved/report_well_imputed.py \\
            --info-files ${inWell_imputed} \\
            --datasets ${datasets} \\
            --output-prefix ${out_prefix} \\
            --rsq-threshold ${params.impute_info_cutoff} \\
            --maf-threshold ${params.maf_thresh}
        """
}

// Plot performance target - chromosome level
process plot_performance_target_chr {
    tag "plot_performance_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/plots_chr/chr${chr}/${ref_panels}", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), file(well_imputed_report), file(well_imputed_report_summary), val(group)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(plot_by_maf)
    script:
        plot_by_maf = "${well_imputed_report.baseName}_chr${chr}.pdf"
        report = well_imputed_report
        xlab = "MAF bins"
        ylab = "Number of well imputed SNPs (Chr ${chr})"
        template "plot_results_by_maf.py"
}

// Duplicate for dataset grouping
process plot_performance_target_chr2 {
    tag "plot_performance_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/plots_chr/chr${chr}/${target_name}", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), file(well_imputed_report), file(well_imputed_report_summary), val(group)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(plot_by_maf)
    script:
        plot_by_maf = "${well_imputed_report.baseName}_chr${chr}.pdf"
        report = well_imputed_report
        xlab = "MAF bins"
        ylab = "Number of well imputed SNPs (Chr ${chr})"
        template "plot_results_by_maf.py"
}

// Report accuracy target - chromosome level
process report_accuracy_target_chr {
    tag "report_acc_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/reports_chr/chr${chr}/${ref_panels}/", overwrite: true, mode:'copy'
    label "medium"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), file(inSNP_acc), val(group)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(outSNP_acc), val(group)
    script:
        outSNP_acc = "${inSNP_acc.baseName}.chr${chr}.imputed_info_report_accuracy"
        datasets = ref_panels.split(',').join(',')
        """
        python3 ${projectDir}/templates/improved/report_accuracy_by_maf.py \\
            --info-files ${inSNP_acc} \\
            --datasets ${datasets} \\
            --output-prefix ${outSNP_acc} \\
            --rsq-threshold ${params.impute_info_cutoff} \\
            --summary-table
        """
}

// Duplicate for dataset grouping
process report_accuracy_target_chr2 {
    tag "report_acc_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/reports_chr/chr${chr}/${target_name}/", overwrite: true, mode:'copy'
    label "medium"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), file(inSNP_acc), val(group)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(outSNP_acc), val(group)
    script:
        outSNP_acc = "${inSNP_acc.baseName}.chr${chr}.imputed_info_report_accuracy"
        datasets = ref_panels.split(',').join(',')
        """
        python3 ${projectDir}/templates/improved/report_accuracy_by_maf.py \\
            --info-files ${inSNP_acc} \\
            --datasets ${datasets} \\
            --output-prefix ${outSNP_acc} \\
            --rsq-threshold ${params.impute_info_cutoff} \\
            --summary-table
        """
}

// Plot accuracy target - chromosome level
process plot_accuracy_target_chr {
    tag "plot_accuracy_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/plots_chr/chr${chr}/${ref_panels}", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), file(accuracy_report), val(group)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(plot_by_maf)
    script:
        plot_by_maf = "${accuracy_report.baseName}_accuracy_by_maf_chr${chr}.pdf"
        report = accuracy_report
        xlab = "MAF bins"
        ylab = "Concordance rate (Chr ${chr})"
        template "plot_results_by_maf.py"
}

// Duplicate for dataset grouping
process plot_accuracy_target_chr2 {
    tag "plot_accuracy_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/plots_chr/chr${chr}/${target_name}", overwrite: true, mode:'copy'
    label "python_plotting"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), file(accuracy_report), val(group)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(plot_by_maf)
    script:
        plot_by_maf = "${accuracy_report.baseName}_accuracy_by_maf_chr${chr}.pdf"
        report = accuracy_report
        xlab = "MAF bins"
        ylab = "Concordance rate (Chr ${chr})"
        template "plot_results_by_maf.py"
}

// Plot R2 vs SNP count - chromosome level
process plot_r2_SNPcount_chr {
    tag "plot_r2_SNPcount_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/plots_chr/chr${chr}/${target_name}", overwrite: true, mode:'copy'
    label "medium"
    label "python_plotting"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), val(infos)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(plot_out)
    script:
        plot_out = "${target_name}_${ref_panels}_chr${chr}_r2_SNPcount.pdf"
        impute_info_cutoff = params.impute_info_cutoff
        template "r2_Frequency_plot.py"
}

// Plot histogram R2 vs SNP count - chromosome level
process plot_hist_r2_SNPcount_chr {
    tag "plot_hist_r2_SNPcount_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/plots_chr/chr${chr}/${target_name}/", overwrite: true, mode:'copy'
    label "medium"
    label "python_plotting"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), val(infos)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(plot_out)
    script:
        plot_out = "${target_name}_${ref_panels}_chr${chr}_r2_SNPcount_hist.pdf"
        impute_info_cutoff = params.impute_info_cutoff
        template "r2_Frequency_plot_histogram.py"
}

// Plot MAF vs R2 - chromosome level
process plot_MAF_r2_chr {
    tag "plot_MAF_r2_${target_name}_${ref_panels}_chr${chr}"
    publishDir "${params.outDir}/plots_chr/chr${chr}/${target_name}", overwrite: true, mode:'copy'
    label "medium"
    label "python_plotting"
    
    input:
        tuple val(target_name), val(ref_panels), val(chr), val(infos)
    output:
        tuple val(target_name), val(ref_panels), val(chr), file(plot_out)
    script:
        plot_out = "${target_name}_${ref_panels}_chr${chr}_MAF_r2.pdf"
        impute_info_cutoff = params.impute_info_cutoff
        template "Frequency_r2_MAF_plot.py"
}