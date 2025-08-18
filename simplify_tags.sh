#!/bin/bash

# Script to simplify verbose Nextflow process tags

# Fix report.nf
sed -i 's/tag "filter_${dataset_name}_${tagName}_${ref_panels.join.*}"/tag "${tagName}"/g' modules/report.nf
sed -i 's/tag "site_by_maf_${dataset_name}"/tag "${dataset_name}"/g' modules/report.nf
sed -i 's/tag "report_wellImputed_${target_name}_${ref_panels.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "plot_performance_dataset_${target_name}_${ref_panels.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "report_acc_${target_name}_${ref_panels.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "plot_accuracy_dataset_${target_name}_${ref_panels.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "frq_${target_name}_${ref_name}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "plot_r2_SNPpos_${target_name}_${ref_name.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "plot_freq_comparison_${target_name}_${ref_name.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "plot_r2_SNPcount_${target_name}_${ref_panels.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "plot_hist_r2_SNPcount_${target_name}_${ref_panels.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "plot_MAF_r2_${target_name}_${ref_panels.*}"/tag "${target_name}"/g' modules/report.nf
sed -i 's/tag "average_r2_${target_name}_${ref_panels.*}"/tag "${target_name}"/g' modules/report.nf

# Fix qc.nf - for chunks, use the full chunk identifier
sed -i 's/tag "${target_name}_${chrm}:${chunk_start}-${chunk_end}"/tag "${target_name}_chr${chrm}_${chunk_start}_${chunk_end}"/g' modules/qc.nf
sed -i 's/tag "${target_name}_${chrm}"/tag "${target_name}"/g' modules/qc.nf

# Fix report_chr.nf
sed -i 's/tag "filter_${dataset_name}_${tagName}_${ref_panels.*}_${chr}"/tag "${tagName}"/g' modules/report_chr.nf
sed -i 's/tag "filter_${dataset_name}_${ref_panels.*}_chr${chr}"/tag "${dataset_name}_chr${chr}"/g' modules/report_chr.nf
sed -i 's/tag "report_wellImputed_${target_name}_${ref_panels.*}_chr${chr}"/tag "${target_name}_chr${chr}"/g' modules/report_chr.nf
sed -i 's/tag "plot_performance_${target_name}_${ref_panels.*}_chr${chr}"/tag "${target_name}_chr${chr}"/g' modules/report_chr.nf
sed -i 's/tag "report_acc_${target_name}_${ref_panels.*}_chr${chr}"/tag "${target_name}_chr${chr}"/g' modules/report_chr.nf
sed -i 's/tag "plot_accuracy_${target_name}_${ref_panels.*}_chr${chr}"/tag "${target_name}_chr${chr}"/g' modules/report_chr.nf
sed -i 's/tag "plot_r2_SNPcount_${target_name}_${ref_panels.*}_chr${chr}"/tag "${target_name}_chr${chr}"/g' modules/report_chr.nf
sed -i 's/tag "plot_hist_r2_SNPcount_${target_name}_${ref_panels.*}_chr${chr}"/tag "${target_name}_chr${chr}"/g' modules/report_chr.nf
sed -i 's/tag "plot_MAF_r2_${target_name}_${ref_panels.*}_chr${chr}"/tag "${target_name}_chr${chr}"/g' modules/report_chr.nf

# Fix subset_vcf.nf
sed -i 's/tag "extract_site_${target_name}_${site_name}"/tag "${target_name}"/g' modules/subset_vcf.nf

echo "Tag simplification complete!"