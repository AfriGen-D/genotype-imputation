#!/usr/bin/env nextflow
nextflow.enable.dsl=2  // Updated to stable DSL2 syntax (from deprecated preview mode)

include { get_chromosome; fill_tags_vcf; check_chromosome; check_files; check_chromosome_vcf; check_mismatch; no_mismatch ; qc_dupl; split_multi_allelic; filter_min_ac; target_qc as target_qc; target_qc as target_qc1; qc_site_missingness as qc_site_missingness1; qc_site_missingness as qc_site_missingness2; sites_only ; combine_vcfs ; combine_infos; combine_csvs as combine_freqs; combine_vcfs_chrm; } from './modules/qc' 
include { vcf_map_simple; extract_site_from_vcf; generate_chunks_vcf; split_target_to_chunk; vcf_map; vcf_freq; info_freq; fill_tags_VCF; sort_vcf ; get_vcf_sites; extract_pop } from './modules/subset_vcf'
include { minimac4_phasing_eagle } from './modules/phasing'
include { impute_minimac4; extract_impute_info; impute_minimac4_1; combineImpute; combineInfo; filter_info_by_target } from './modules/impute'
include { filter_info; report_site_by_maf; plot_freq_comparison; report_well_imputed_by_target; plot_performance_target; 
report_accuracy_target; plot_accuracy_target; generate_frequency; plot_r2_SNPpos; plot_r2_SNPcount; plot_hist_r2_SNPcount; plot_MAF_r2; 
average_r2 } from './modules/report'
include { filter_info_by_target_chr; filter_info_by_target_chr2; report_well_imputed_by_target_chr; report_well_imputed_by_target_chr2;
plot_performance_target_chr; plot_performance_target_chr2; report_accuracy_target_chr; report_accuracy_target_chr2;
plot_accuracy_target_chr; plot_accuracy_target_chr2; plot_r2_SNPcount_chr; plot_hist_r2_SNPcount_chr; plot_MAF_r2_chr } from './modules/report_chr'
include { run_qc_plots } from './modules/qc_plots'


// Header log info
def intro(){
log.info ""
log.info """
=======================================================
h3achipimputation v${params.version}"
======================================================= """
    def summary = [:]
    summary['Pipeline Name']    = 'h3achipimputation'
    summary['Pipeline version'] = params.version
    summary['Run Name']         = workflow.runName
    summary['Target datasets']  = params.target_datasets.collect{ "${it[0]} (${it[1]})" }.join(', ')
    summary['Reference panels'] = params.ref_panels.collect{ "${it[0]} (${it[2]})" }.join(', ')
    summary['Max Memory']       = params.max_memory
    summary['Max CPUs']         = params.max_cpus
    summary['Max Time']         = params.max_time
    summary['Output dir']       = params.outDir
    summary['Working dir']      = workflow.workDir
    summary['Script dir']       = workflow.projectDir
    summary['Current path']     = "$PWD"
    summary['Git info']         = "${workflow.repository} - ${workflow.revision} [${workflow.commitId}]"
    summary['Command line']     = workflow.commandLine
    if(workflow.containerEngine) {
        summary['Container Engine'] = workflow.containerEngine
        summary['Container'] = workflow.container
        summary['Current home'] = "$HOME"
        summary['Current user'] = "$USER"
        summary['Current path'] = "$PWD"
        summary['Working dir'] = workflow.workDir
        summary['Output dir'] = params.outDir
        summary['Script dir'] = workflow.projectDir
        summary['Config Profile'] = workflow.profile
    }

    if(params.email) summary['E-mail Address'] = params.email
    log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
    log.info "======================================================="
    log.info ""
}

workflow preprocess {
    take: datasets

    main:
        //// check if study genotype files exist
        target_datasets = []
        datasets.each { name, vcf ->
            check_files([vcf])
            target_datasets << [name, file(vcf)]
        }
        target_datasets = Channel.from(target_datasets)

        //// check if eagle map file exists
        if(params.eagle_genetic_map) {
            check_files([params.eagle_genetic_map])
        }

        //// check if fasta reference genome and index
        if(params.reference_genome) {
            check_files([params.reference_genome, "${params.reference_genome}.fai"])
        }

        // //// Check chromosome
        // get_chromosome( target_datasets )
        // chromosomes = get_chromosome.out
        //     .map{ dataset, dataset_vcf, map_file -> check_chromosome_vcf(dataset, dataset_vcf, map_file, params.chromosomes) }
        //     .map{ dataset, dataset_vcf, map_file, chrms -> chrms.unique() }
        
        //// check if reference panel files exist - handle ALL chromosomes
        def resolved_chromosomes = []
        if (params.chromosomes == 'ALL' || params.chromosomes == '') {
            // When ALL is specified, use both b37 and b38 chromosome formats
            // This prevents the chrALL file path error while allowing the pipeline to continue
            log.info "|-- INFO: Chromosomes set to 'ALL' - checking reference files for common chromosomes (both b37 and b38 formats)"
            
            // Try both b37 (1, 2, 3...) and b38 (chr1, chr2, chr3...) formats
            def b37_chrms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
            def b38_chrms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
            
            // Check which format the reference panel files use by testing a few chromosomes
            def test_chrm_found = false
            params.ref_panels.each{ ref_name, ref_m3vcf, ref_vcf ->
                if (!test_chrm_found) {
                    // Test chr1 first (b38 format)
                    def test_vcf_b38 = sprintf(ref_vcf, 'chr1')
                    def test_m3vcf_b38 = sprintf(ref_m3vcf, 'chr1')
                    if (file(test_vcf_b38).exists() || file(test_m3vcf_b38).exists()) {
                        resolved_chromosomes = b38_chrms
                        test_chrm_found = true
                        log.info "|-- INFO: Detected b38 chromosome format (chr1, chr2, ...) in reference panels"
                    } else {
                        // Test 1 (b37 format)
                        def test_vcf_b37 = sprintf(ref_vcf, '1')
                        def test_m3vcf_b37 = sprintf(ref_m3vcf, '1')
                        if (file(test_vcf_b37).exists() || file(test_m3vcf_b37).exists()) {
                            resolved_chromosomes = b37_chrms
                            test_chrm_found = true
                            log.info "|-- INFO: Detected b37 chromosome format (1, 2, ...) in reference panels"
                        }
                    }
                }
            }
            
            // If no files found, default to both formats for checking
            if (!test_chrm_found) {
                resolved_chromosomes = b37_chrms + b38_chrms
                log.warn "|-- WARN: Could not determine chromosome format from reference panels, will check both b37 and b38 formats"
            }
        } else {
            resolved_chromosomes = params.chromosomes.split(',')
        }
        
        // Only check files for chromosomes that have corresponding reference files
        resolved_chromosomes.each{ chrm ->
            params.ref_panels.each{ ref_name, ref_m3vcf, ref_vcf ->
                vcf = sprintf(ref_vcf, chrm)
                m3vcf = sprintf(ref_m3vcf, chrm)
                
                // Only check files if they should exist (skip if file patterns don't make sense)
                if (!vcf.contains('chrALL') && !m3vcf.contains('chrALL')) {
                    if(vcf.endsWith("vcf.gz")){
                        vcf_idx = "${vcf}.tbi" 
                    }
                    else if(vcf.endsWith("bcf")){
                        vcf_idx = "${vcf}.csi" 
                    }
                    try {
                        check_files([ m3vcf, vcf, vcf_idx ])
                    } catch (Exception e) {
                        log.warn "|-- WARN: Reference files for chromosome ${chrm} not found: ${m3vcf}, ${vcf}"
                        log.warn "|--       Pipeline will determine actual available chromosomes from target data"
                    }
                } else {
                    log.warn "|-- WARN: Skipping file check for chromosome '${chrm}' as it resolves to invalid file path"
                }
            }
        }

        ////// QC
        check_mismatch(target_datasets.map{ dataset, dataset_vcf -> [ dataset, '', '', '', file(dataset_vcf), file(params.reference_genome) ] } )
        qc_dupl( target_datasets.map{ dataset, dataset_vcf -> [ dataset, '', '', '', file(dataset_vcf) ] } )
        split_multi_allelic(qc_dupl.out)
        fill_tags_vcf(split_multi_allelic.out)
        filter_min_ac(fill_tags_vcf.out.map{ dataset, chrm, start, end, vcf -> [ dataset, chrm, start, end, file(vcf), " --min-ac ${params.min_ac} --max-alleles ${params.max_alleles} --min-alleles ${params.min_alleles} -v snps "  ] })

    emit:
        // chrms = chromosomes
        dataset_qc = filter_min_ac.out
        // get_chromosome = get_chromosome.out
}

workflow subset{
    take: data
    
    main:
        get_chromosome( data.map{ dataset, chrm, start, end, vcf, map_file -> [ dataset, file(vcf) ] } )
        data_chrms = get_chromosome.out
            .map{ dataset, dataset_vcf, map_file -> check_chromosome_vcf(dataset, dataset_vcf, map_file, params.chromosomes) }
            .map{ dataset, dataset_vcf, map_file, chrms -> [ dataset, file(dataset_vcf), file(map_file), chrms.unique().join(',') ] }
        
        generate_chunks_vcf(data_chrms.map{ dataset, vcf, map_file, chrms -> [ dataset, file(vcf), file(map_file), chrms, params.chunk_size ] })
        chunks_datas = generate_chunks_vcf.out.flatMap{ dataset, vcf, chunk_file ->
            datas = []
            chunks = file(chunk_file).readLines()
            chunks.each{ chunk_data ->
                data = chunk_data.trim().split(',')
                chrm = data[0]
                chunk_start = data[1]
                chunk_end = data[2]
                datas << [dataset, chrm, chunk_start, chunk_end, dataset, file(vcf)]
            }
            return datas
        }
        split_target_to_chunk(chunks_datas)

    emit: 
        chunks = split_target_to_chunk.out
}

workflow phasing{
    take: data
    
    main:
        minimac4_phasing_eagle(data)

    emit: 
        chunks_phased = minimac4_phasing_eagle.out
}

workflow impute{
    take: data
    
    main:
        impute_minimac4(data)
        extract_impute_info(impute_minimac4.out)

    emit: 
        data
        chunks_imputed = extract_impute_info.out
}

workflow report_by_ref{
    take: data

    main:
        /// /// By Refecene panel
        // combine by chrom, dataset, refpanel
        imputeCombine_ref = data
                .groupTuple( by:[1] )
                .map{ datasets, refpanel, vcfs, imputed_vcfs, imputed_infos -> [ refpanel, datasets.join(','), '', imputed_infos.join(',') ] }
        filter_info_by_target( imputeCombine_ref )

        /// change to group_by_maf
        report_well_imputed_by_target( filter_info_by_target.out.map{ target_name, ref_panels, wellInfo, accInfo -> [ target_name, ref_panels, file(wellInfo) ]} )
        
        //// Plot performance all targets by maf for a reference panel
        plot_performance_target( report_well_imputed_by_target.out.map{ target_name, ref_panels, wellInfo, wellInfo_summary -> [ target_name, ref_panels, file(wellInfo), file(wellInfo_summary), 'DATASETS' ]} )

        //// Accuracy/Concordance
        report_accuracy_target( filter_info_by_target.out.map{ target_name, ref_panels, wellInfo, accInfo -> [ target_name, ref_panels, file(accInfo), 'DATASETS' ]} )
        plot_accuracy_target ( report_accuracy_target.out )
    emit:
        data
}

workflow report_by_dataset{
    take: data

    main:
        /// /// By Dataset
        // combine by chrom, dataset, refpanel
        
        imputeCombine_ref = data
                .groupTuple( by:[0] )
                .map{ dataset, refpanels, vcfs, imputed_vcfs, imputed_infos -> [ dataset, refpanels.join(','), '', imputed_infos.join(',') ] }
        filter_info_by_target( imputeCombine_ref )

        ///// Number of well imputed snps
        /// change to group_by_maf
        report_well_imputed_by_target( filter_info_by_target.out.map{ target_name, ref_panels, wellInfo, accInfo -> [ target_name, ref_panels, file(wellInfo) ]} )
        /// Plot performance all targets by maf for a reference panel
        plot_performance_target( report_well_imputed_by_target.out.map{ target_name, ref_panels, wellInfo, wellInfo_summary -> [ target_name, ref_panels, file(wellInfo), file(wellInfo_summary), 'REFERENCE_PANELS' ]} )

        //// Accuracy/Concordance
        report_accuracy_target( filter_info_by_target.out.map{ target_name, ref_panels, wellInfo, accInfo -> [ target_name, ref_panels, file(accInfo), 'REFERENCE_PANELS' ]} )
        plot_accuracy_target ( report_accuracy_target.out )

        // Plot number of imputed SNPs over the mean r2 for all reference panels
        input = imputeCombine_ref
        .map{ dataset, refpanels, chrm, infos -> [dataset, refpanels, infos]}

        // Plot number of imputed SNPs over the mean r2 for all reference panels
        plot_r2_SNPcount(input)

        // Plot histograms of number of imputed SNPs over the mean r2 for all reference panels
        plot_hist_r2_SNPcount(input)

        // Plot MAF of imputed SNPs over r2 for all references
        plot_MAF_r2(input)

    emit:
        data
}

workflow report_by_ref_chromosome {
    take: data_with_chr

    main:
        /// By Reference panel - chromosome level
        // Group by chromosome and reference panel
        imputeCombine_ref_chr = data_with_chr
                .groupTuple( by:[0, 2] )  // Group by chromosome and ref_panel
                .map{ chr, datasets, refpanel, vcfs, imputed_vcfs, imputed_infos -> 
                    [ refpanel, datasets.join(','), chr, imputed_infos.join(',') ] 
                }
        
        filter_info_by_target_chr( imputeCombine_ref_chr )

        // Plot performance by chromosome
        report_well_imputed_by_target_chr( 
            filter_info_by_target_chr.out.map{ target_name, ref_panels, chr, wellInfo, accInfo -> 
                [ target_name, ref_panels, chr, file(wellInfo) ]
            } 
        )
        
        plot_performance_target_chr( 
            report_well_imputed_by_target_chr.out.map{ target_name, ref_panels, chr, wellInfo, wellInfo_summary -> 
                [ target_name, ref_panels, chr, file(wellInfo), file(wellInfo_summary), 'DATASETS' ]
            } 
        )

        // Accuracy/Concordance by chromosome
        report_accuracy_target_chr( 
            filter_info_by_target_chr.out.map{ target_name, ref_panels, chr, wellInfo, accInfo -> 
                [ target_name, ref_panels, chr, file(accInfo), 'DATASETS' ]
            } 
        )
        plot_accuracy_target_chr( report_accuracy_target_chr.out )

    emit:
        data = filter_info_by_target_chr.out
}

workflow report_by_dataset_chromosome {
    take: data_with_chr

    main:
        /// By Dataset - chromosome level
        // Group by chromosome and dataset
        imputeCombine_dataset_chr = data_with_chr
                .groupTuple( by:[0, 1] )  // Group by chromosome and dataset
                .map{ chr, dataset, refpanels, vcfs, imputed_vcfs, imputed_infos -> 
                    [ dataset, refpanels.join(','), chr, imputed_infos.join(',') ] 
                }
        
        filter_info_by_target_chr2( imputeCombine_dataset_chr )

        // Plot performance by chromosome
        report_well_imputed_by_target_chr2( 
            filter_info_by_target_chr2.out.map{ target_name, ref_panels, chr, wellInfo, accInfo -> 
                [ target_name, ref_panels, chr, file(wellInfo) ]
            } 
        )
        
        plot_performance_target_chr2( 
            report_well_imputed_by_target_chr2.out.map{ target_name, ref_panels, chr, wellInfo, wellInfo_summary -> 
                [ target_name, ref_panels, chr, file(wellInfo), file(wellInfo_summary), 'REFERENCE_PANELS' ]
            } 
        )

        // Accuracy/Concordance by chromosome
        report_accuracy_target_chr2( 
            filter_info_by_target_chr2.out.map{ target_name, ref_panels, chr, wellInfo, accInfo -> 
                [ target_name, ref_panels, chr, file(accInfo), 'REFERENCE_PANELS' ]
            } 
        )
        plot_accuracy_target_chr2( report_accuracy_target_chr2.out )

        // R2 and MAF plots by chromosome
        plot_r2_SNPcount_chr( imputeCombine_dataset_chr.map{ dataset, refpanels, chr, infos -> 
            [dataset, refpanels, chr, infos]
        })
        
        plot_hist_r2_SNPcount_chr( imputeCombine_dataset_chr.map{ dataset, refpanels, chr, infos -> 
            [dataset, refpanels, chr, infos]
        })
        
        plot_MAF_r2_chr( imputeCombine_dataset_chr.map{ dataset, refpanels, chr, infos -> 
            [dataset, refpanels, chr, infos]
        })

    emit:
        data = filter_info_by_target_chr2.out
}

workflow {

    intro()

    //// Data preparation
    preprocess(params.target_datasets)
    
    // //// Data chunking
    subset(preprocess.out.dataset_qc)

    //// Phasing
    phasing_data = Channel.from(params.ref_panels)
        .combine(subset.out.chunks)
        .flatMap{ ref_name, ref_m3vcf, ref_vcf, dataset, chrm, start, end, dataset_, dataset_vcf ->
            vcf = sprintf(ref_vcf, chrm)
            m3vcf = sprintf(ref_m3vcf, chrm)
            if(vcf.endsWith("vcf.gz")){
                vcf_idx = "${vcf}.tbi" 
            }
            else if(vcf.endsWith("bcf")){
                vcf_idx = "${vcf}.csi" 
            }
            return [ [ chrm, ref_name, file(m3vcf), file(vcf), file(vcf_idx), file(params.eagle_genetic_map), start, end, dataset, dataset, file(dataset_vcf)] ]
        }
    phasing(phasing_data)

    //// Imputation
    // phasing.out.chunks_phased.view()
    impute(phasing.out.chunks_phased)
    
    // Reporting
    // impute.out.chunks_imputed.view()
    impute_data = impute.out.chunks_imputed
                .map{chr, fwd, rev, test_data, ref, imputed_vcf, 
                imputed_info, tst_data -> [test_data, ref, imputed_vcf, imputed_info]}
                .combine(params.target_datasets, by:0)
                .map {test_data, ref, imputed_vcf, imputed_info, orig_vcf 
                -> [test_data, ref, orig_vcf, imputed_vcf, imputed_info]}

    // //// GENOME-WIDE REPORTING
    // //// Report by Reference - Genome-wide
    report_by_ref( impute_data )

    // //// Report by datasets - Genome-wide
    report_by_dataset( impute_data )
    
    // //// CHROMOSOME-LEVEL REPORTING
    // Prepare data with chromosome information
    impute_data_chr_ref = impute.out.chunks_imputed
        .map{chr, fwd, rev, test_data, ref, imputed_vcf, imputed_info, tst_data -> 
            [chr, test_data, ref, imputed_vcf, imputed_info]}
        .combine(params.target_datasets, by:1)
        .map {test_data, chr, ref, imputed_vcf, imputed_info, orig_vcf -> 
            [chr, test_data, ref, orig_vcf, imputed_vcf, imputed_info]}
    
    // //// Report by Reference - Chromosome level
    report_by_ref_chromosome( impute_data_chr_ref )
    
    // //// Report by datasets - Chromosome level  
    report_by_dataset_chromosome( impute_data_chr_ref )

    // // Generate dataset frequencies - need chromosome for ref panel path
    impute_data_with_chr = impute.out.chunks_imputed
                .map{chr, fwd, rev, test_data, ref, imputed_vcf, 
                imputed_info, tst_data -> [chr, test_data, ref, imputed_vcf, imputed_info]}
    
    input = impute_data_with_chr
    .map{ chr, target_name, ref_name, impute_vcf, info -> 
        // Get the reference panel VCF path and replace %s with chromosome
        def ref_panel = params.ref_panels.find { it[0] == ref_name }
        if (ref_panel) {
            def ref_vcf_path = ref_panel[2].replace('%s', chr)
            [target_name, ref_name, file(impute_vcf), file(ref_vcf_path)]
        }
    }
    .filter { it != null }  // Remove any null entries
    generate_frequency(input)

    // // Plot frequency Comparison
    freq_comp = impute_data_with_chr.map {chr, target_name, ref_name, impute_vcf, info -> 
    [target_name, ref_name, info]}
    .combine(generate_frequency.out, by:[0,1])
    plot_freq_comparison(freq_comp)

    // // Plot number of imputed SNPs over the mean r2 for all reference panels
    combineInfo_frq = impute_data_with_chr.map{ chr, target_name, ref_name, impute_vcf, info ->[ target_name, ref_name, info, params.maf_thresh]}
    .combine(generate_frequency.out, by:[0,1])
    .map { target_name, ref_name, info, maf_thresh, target_frq, ref_frq -> 
    [target_name, ref_name, info, maf_thresh, target_frq]}
    plot_r2_SNPpos(combineInfo_frq)

    // // compute for average rsquared values
    rsquared_input = impute_data.map{ target_name, ref_name, vcf, impute_vcf, info ->[ target_name, ref_name, info]}
    average_r2(rsquared_input)
}