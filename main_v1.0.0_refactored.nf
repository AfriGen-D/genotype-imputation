#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================================
// H3A CHIP IMPUTATION PIPELINE v1.0.0
// ================================================================================
// Refactored: 16-aug-2025
// Maintainer: H3ABioNet Informatics
// ================================================================================

// ================================================================================
// MODULE IMPORTS
// ================================================================================

// Quality Control Modules
include { 
    get_chromosome; fill_tags_vcf; check_chromosome; check_files; 
    check_chromosome_vcf; check_mismatch; no_mismatch; qc_dupl; 
    split_multi_allelic; filter_min_ac; 
    target_qc as target_qc; target_qc as target_qc1; 
    qc_site_missingness as qc_site_missingness1; 
    qc_site_missingness as qc_site_missingness2; 
    sites_only; combine_vcfs; combine_infos; 
    combine_csvs as combine_freqs; combine_vcfs_chrm 
} from './modules/qc' 

// VCF Processing Modules
include { 
    vcf_map_simple; extract_site_from_vcf; generate_chunks_vcf; 
    split_target_to_chunk; vcf_map; vcf_freq; info_freq; 
    fill_tags_VCF; sort_vcf; get_vcf_sites; extract_pop 
} from './modules/subset_vcf'

// Phasing Modules
include { minimac4_phasing_eagle } from './modules/phasing'

// Imputation Modules
include { 
    impute_minimac4; extract_impute_info; impute_minimac4_1; 
    combineImpute; combineInfo; filter_info_by_target 
} from './modules/impute'

// Reporting Modules (Genome-wide)
include { 
    filter_info; report_site_by_maf; plot_freq_comparison; 
    report_well_imputed_by_target; plot_performance_target; 
    report_accuracy_target; plot_accuracy_target; generate_frequency; 
    plot_r2_SNPpos; plot_r2_SNPcount; plot_hist_r2_SNPcount; 
    plot_MAF_r2; average_r2 
} from './modules/report'

// Reporting Modules (Chromosome-level)
include { 
    filter_info_by_target_chr; filter_info_by_target_chr2; 
    report_well_imputed_by_target_chr; report_well_imputed_by_target_chr2;
    plot_performance_target_chr; plot_performance_target_chr2; 
    report_accuracy_target_chr; report_accuracy_target_chr2;
    plot_accuracy_target_chr; plot_accuracy_target_chr2; 
    plot_r2_SNPcount_chr; plot_hist_r2_SNPcount_chr; plot_MAF_r2_chr 
} from './modules/report_chr'

// QC Plotting Modules
include { run_qc_plots } from './modules/qc_plots'

// ================================================================================
// HELPER FUNCTIONS
// ================================================================================

/**
 * Display pipeline header and configuration summary
 */
def displayHeader() {
    log.info ""
    log.info """
    ╔════════════════════════════════════════════════════════════╗
    ║                H3A CHIP IMPUTATION PIPELINE                ║
    ║                    Version ${params.version} (16-aug-2025)                ║
    ╚════════════════════════════════════════════════════════════╝
    """.stripIndent()
    
    def summary = [:]
    summary['Pipeline Name']    = 'h3achipimputation'
    summary['Pipeline Version'] = params.version
    summary['Build Date']       = '16-aug-2025'
    summary['Run Name']         = workflow.runName
    summary['Target Datasets']  = params.target_datasets.collect{ "${it[0]} (${it[1]})" }.join(', ')
    summary['Reference Panels'] = params.ref_panels.collect{ "${it[0]} (${it[2]})" }.join(', ')
    summary['Chromosomes']      = params.chromosomes ?: 'ALL'
    summary['Chunk Size']       = params.chunk_size ?: 'Default'
    summary['Min AC Filter']    = params.min_ac
    summary['Max Memory']       = params.max_memory
    summary['Max CPUs']         = params.max_cpus
    summary['Max Time']         = params.max_time
    summary['Output Directory'] = params.outDir
    summary['Working Directory']= workflow.workDir
    summary['Script Directory'] = workflow.projectDir
    summary['Current Path']     = "$PWD"
    
    if(workflow.repository) {
        summary['Git Repository'] = workflow.repository
        summary['Git Revision']   = workflow.revision
        summary['Git Commit ID']  = workflow.commitId
    }
    
    summary['Command Line'] = workflow.commandLine
    
    if(workflow.containerEngine) {
        summary['Container Engine'] = workflow.containerEngine
        summary['Container']        = workflow.container
        summary['Config Profile']   = workflow.profile
    }
    
    if(params.email) {
        summary['Email Address'] = params.email
    }
    
    log.info "Configuration Summary:"
    log.info "─────────────────────────────────────────────────────────────"
    summary.each { k, v -> 
        log.info "${k.padRight(20)}: $v"
    }
    log.info "═════════════════════════════════════════════════════════════"
    log.info ""
}

/**
 * Resolve chromosome naming convention (b37 vs b38)
 * @return List of chromosomes to process
 */
def resolveChromosomes() {
    def resolved_chromosomes = []
    
    if (params.chromosomes == 'ALL' || params.chromosomes == '') {
        log.info "→ Detecting chromosome naming convention..."
        
        // Define both naming conventions
        def b37_chrms = (1..22).collect { it.toString() }
        def b38_chrms = (1..22).collect { "chr${it}" }
        
        // Detect which format is used in reference panels
        def format_detected = detectChromosomeFormat()
        
        if (format_detected == 'b38') {
            resolved_chromosomes = b38_chrms
            log.info "  ✓ Using b38 format (chr1, chr2, ...)"
        } else if (format_detected == 'b37') {
            resolved_chromosomes = b37_chrms
            log.info "  ✓ Using b37 format (1, 2, ...)"
        } else {
            // If detection fails, try both formats
            resolved_chromosomes = b37_chrms + b38_chrms
            log.warn "  ⚠ Could not detect format, will try both b37 and b38"
        }
    } else {
        // Use explicitly specified chromosomes
        resolved_chromosomes = params.chromosomes.split(',')
        log.info "→ Using specified chromosomes: ${resolved_chromosomes.join(', ')}"
    }
    
    return resolved_chromosomes
}

/**
 * Detect chromosome format from reference panel files
 * @return String indicating format ('b37', 'b38', or 'unknown')
 */
def detectChromosomeFormat() {
    for (panel in params.ref_panels) {
        def (ref_name, ref_m3vcf, ref_vcf) = panel
        
        // Test b38 format first
        def test_b38 = sprintf(ref_vcf, 'chr1')
        if (file(test_b38).exists() || file(sprintf(ref_m3vcf, 'chr1')).exists()) {
            return 'b38'
        }
        
        // Test b37 format
        def test_b37 = sprintf(ref_vcf, '1')
        if (file(test_b37).exists() || file(sprintf(ref_m3vcf, '1')).exists()) {
            return 'b37'
        }
    }
    
    return 'unknown'
}

/**
 * Validate reference panel files for specified chromosomes
 * @param chromosomes List of chromosomes to validate
 */
def validateReferenceFiles(chromosomes) {
    log.info "→ Validating reference panel files..."
    
    def valid_count = 0
    def missing_count = 0
    
    chromosomes.each { chrm ->
        params.ref_panels.each { ref_name, ref_m3vcf, ref_vcf ->
            def vcf = sprintf(ref_vcf, chrm)
            def m3vcf = sprintf(ref_m3vcf, chrm)
            
            // Skip invalid chromosome placeholders
            if (vcf.contains('chrALL') || m3vcf.contains('chrALL')) {
                return
            }
            
            // Determine index file type
            def vcf_idx = vcf.endsWith("vcf.gz") ? "${vcf}.tbi" : "${vcf}.csi"
            
            // Check if files exist
            def files_exist = file(m3vcf).exists() && 
                             file(vcf).exists() && 
                             file(vcf_idx).exists()
            
            if (files_exist) {
                valid_count++
            } else {
                missing_count++
                log.debug "  Missing: ${ref_name} chromosome ${chrm}"
            }
        }
    }
    
    log.info "  ✓ Found ${valid_count} valid reference files"
    if (missing_count > 0) {
        log.warn "  ⚠ Missing ${missing_count} reference files (will be skipped)"
    }
}

/**
 * Validate input files and parameters
 */
def validateInputs() {
    log.info "→ Validating input files and parameters..."
    
    // Check target datasets
    params.target_datasets.each { name, vcf ->
        if (!file(vcf).exists()) {
            error "Target dataset file not found: ${vcf}"
        }
    }
    
    // Check eagle genetic map
    if (params.eagle_genetic_map && !file(params.eagle_genetic_map).exists()) {
        error "Eagle genetic map file not found: ${params.eagle_genetic_map}"
    }
    
    // Check reference genome
    if (params.reference_genome) {
        if (!file(params.reference_genome).exists()) {
            error "Reference genome file not found: ${params.reference_genome}"
        }
        if (!file("${params.reference_genome}.fai").exists()) {
            error "Reference genome index not found: ${params.reference_genome}.fai"
        }
    }
    
    log.info "  ✓ All required input files validated"
}

// ================================================================================
// SUB-WORKFLOWS
// ================================================================================

/**
 * Data preprocessing workflow
 * - Quality control
 * - Variant filtering
 * - File validation
 */
workflow preprocess {
    take: 
        datasets  // Input dataset specifications
    
    main:
        // Validate all input files
        validateInputs()
        
        // Create channel from target datasets
        target_datasets = Channel.from(datasets.collect { name, vcf -> 
            [name, file(vcf)] 
        })
        
        // Resolve and validate chromosomes
        def chromosomes = resolveChromosomes()
        validateReferenceFiles(chromosomes)
        
        // Quality control pipeline
        check_mismatch(
            target_datasets.map{ dataset, vcf -> 
                [dataset, '', '', '', file(vcf), file(params.reference_genome)]
            }
        )
        
        qc_dupl(
            target_datasets.map{ dataset, vcf -> 
                [dataset, '', '', '', file(vcf)]
            }
        )
        
        split_multi_allelic(qc_dupl.out)
        fill_tags_vcf(split_multi_allelic.out)
        
        filter_min_ac(
            fill_tags_vcf.out.map{ dataset, chrm, start, end, vcf -> 
                [dataset, chrm, start, end, file(vcf), 
                 "--min-ac ${params.min_ac} --max-alleles ${params.max_alleles} --min-alleles ${params.min_alleles} -v snps"]
            }
        )
    
    emit:
        dataset_qc = filter_min_ac.out
}

/**
 * Data subsetting workflow
 * - Chromosome extraction
 * - Chunk generation
 * - Target splitting
 */
workflow subset {
    take: 
        data  // QC'd dataset
    
    main:
        // Extract chromosomes from VCF
        get_chromosome(
            data.map{ dataset, chrm, start, end, vcf, map_file -> 
                [dataset, file(vcf)]
            }
        )
        
        // Process chromosome information
        data_chrms = get_chromosome.out
            .map{ dataset, vcf, map_file -> 
                check_chromosome_vcf(dataset, vcf, map_file, params.chromosomes)
            }
            .map{ dataset, vcf, map_file, chrms -> 
                [dataset, file(vcf), file(map_file), chrms.unique().join(',')]
            }
        
        // Generate chunks for parallel processing
        generate_chunks_vcf(
            data_chrms.map{ dataset, vcf, map_file, chrms -> 
                [dataset, file(vcf), file(map_file), chrms, params.chunk_size]
            }
        )
        
        // Split data into chunks
        chunks_datas = generate_chunks_vcf.out.flatMap{ dataset, vcf, chunk_file ->
            file(chunk_file).readLines().collect{ chunk_data ->
                def (chrm, start, end) = chunk_data.trim().split(',')
                [dataset, chrm, start, end, dataset, file(vcf)]
            }
        }
        
        split_target_to_chunk(chunks_datas)
    
    emit: 
        chunks = split_target_to_chunk.out
}

/**
 * Phasing workflow using Eagle
 */
workflow phasing {
    take: 
        data  // Chunked data with reference panel info
    
    main:
        minimac4_phasing_eagle(data)
    
    emit: 
        chunks_phased = minimac4_phasing_eagle.out
}

/**
 * Imputation workflow using Minimac4
 */
workflow impute {
    take: 
        data  // Phased data
    
    main:
        impute_minimac4(data)
        extract_impute_info(impute_minimac4.out)
    
    emit: 
        chunks_imputed = extract_impute_info.out
}

/**
 * Generate comprehensive reports by reference panel
 */
workflow report_by_ref {
    take: 
        data  // Imputed data
    
    main:
        // Group data by reference panel
        grouped_data = data
            .groupTuple(by: [1])
            .map{ datasets, refpanel, vcfs, imputed_vcfs, imputed_infos -> 
                def num_files = imputed_infos.size()
                def datasets_replicated = ([refpanel] * num_files).join(',')
                [refpanel, datasets_replicated, '', imputed_infos.join(',')]
            }
        
        filter_info_by_target(grouped_data)
        
        // Generate well-imputed SNPs report
        report_well_imputed_by_target(
            filter_info_by_target.out.map{ target, refs, wellInfo, accInfo -> 
                [target, refs, file(wellInfo)]
            }
        )
        
        // Generate performance plots
        plot_performance_target(
            report_well_imputed_by_target.out.map{ target, refs, wellInfo, summary -> 
                [target, refs, file(wellInfo), file(summary), 'DATASETS']
            }
        )
        
        // Generate accuracy reports
        report_accuracy_target(
            filter_info_by_target.out.map{ target, refs, wellInfo, accInfo -> 
                [target, refs, file(accInfo), 'DATASETS']
            }
        )
        
        plot_accuracy_target(report_accuracy_target.out)
    
    emit:
        reports = filter_info_by_target.out
}

/**
 * Generate comprehensive reports by dataset
 */
workflow report_by_dataset {
    take: 
        data  // Imputed data
    
    main:
        // Group data by dataset
        grouped_data = data
            .groupTuple(by: [0])
            .map{ dataset, refpanels, vcfs, imputed_vcfs, imputed_infos -> 
                def num_files = imputed_infos.size()
                def datasets_replicated = ([dataset] * num_files).join(',')
                [dataset, datasets_replicated, '', imputed_infos.join(',')]
            }
        
        filter_info_by_target(grouped_data)
        
        // Generate well-imputed SNPs report
        report_well_imputed_by_target(
            filter_info_by_target.out.map{ target, refs, wellInfo, accInfo -> 
                [target, refs, file(wellInfo)]
            }
        )
        
        // Generate performance plots
        plot_performance_target(
            report_well_imputed_by_target.out.map{ target, refs, wellInfo, summary -> 
                [target, refs, file(wellInfo), file(summary), 'REFERENCE_PANELS']
            }
        )
        
        // Generate accuracy reports
        report_accuracy_target(
            filter_info_by_target.out.map{ target, refs, wellInfo, accInfo -> 
                [target, refs, file(accInfo), 'REFERENCE_PANELS']
            }
        )
        
        plot_accuracy_target(report_accuracy_target.out)
        
        // Additional analysis plots
        input_data = grouped_data.map{ dataset, refpanels, chrm, infos -> 
            [dataset, refpanels, infos]
        }
        
        plot_r2_SNPcount(input_data)
        plot_hist_r2_SNPcount(input_data)
        plot_MAF_r2(input_data)
    
    emit:
        reports = filter_info_by_target.out
}

/**
 * Generate chromosome-level reports by reference panel
 */
workflow report_by_ref_chromosome {
    take: 
        data_with_chr  // Imputed data with chromosome info
    
    main:
        // Process chromosome-level data
        chr_data = data_with_chr
            .map{ chr, dataset, refpanel, vcf, imputed_vcf, imputed_info ->
                [dataset, refpanel, chr, imputed_info.toString()]
            }
        
        filter_info_by_target_chr(chr_data)
        
        // Generate reports
        report_well_imputed_by_target_chr(
            filter_info_by_target_chr.out.map{ target, refs, chr, wellInfo, accInfo -> 
                [target, refs, chr, file(wellInfo)]
            }
        )
        
        plot_performance_target_chr(
            report_well_imputed_by_target_chr.out.map{ target, refs, chr, wellInfo, summary -> 
                [target, refs, chr, file(wellInfo), file(summary), 'DATASETS']
            }
        )
        
        report_accuracy_target_chr(
            filter_info_by_target_chr.out.map{ target, refs, chr, wellInfo, accInfo -> 
                [target, refs, chr, file(accInfo), 'DATASETS']
            }
        )
        
        plot_accuracy_target_chr(report_accuracy_target_chr.out)
    
    emit:
        reports = filter_info_by_target_chr.out
}

/**
 * Generate chromosome-level reports by dataset
 */
workflow report_by_dataset_chromosome {
    take: 
        data_with_chr  // Imputed data with chromosome info
    
    main:
        // Process chromosome-level data
        chr_data = data_with_chr
            .map{ chr, dataset, refpanel, vcf, imputed_vcf, imputed_info ->
                [dataset, refpanel, chr, imputed_info.toString()]
            }
        
        filter_info_by_target_chr2(chr_data)
        
        // Generate reports
        report_well_imputed_by_target_chr2(
            filter_info_by_target_chr2.out.map{ target, refs, chr, wellInfo, accInfo -> 
                [target, refs, chr, file(wellInfo)]
            }
        )
        
        plot_performance_target_chr2(
            report_well_imputed_by_target_chr2.out.map{ target, refs, chr, wellInfo, summary -> 
                [target, refs, chr, file(wellInfo), file(summary), 'REFERENCE_PANELS']
            }
        )
        
        report_accuracy_target_chr2(
            filter_info_by_target_chr2.out.map{ target, refs, chr, wellInfo, accInfo -> 
                [target, refs, chr, file(accInfo), 'REFERENCE_PANELS']
            }
        )
        
        plot_accuracy_target_chr2(report_accuracy_target_chr2.out)
        
        // Additional chromosome-level plots
        plot_r2_SNPcount_chr(
            chr_data.map{ dataset, refpanels, chr, infos -> 
                [dataset, refpanels, chr, infos]
            }
        )
        
        plot_hist_r2_SNPcount_chr(
            chr_data.map{ dataset, refpanels, chr, infos -> 
                [dataset, refpanels, chr, infos]
            }
        )
        
        plot_MAF_r2_chr(
            chr_data.map{ dataset, refpanels, chr, infos -> 
                [dataset, refpanels, chr, infos]
            }
        )
    
    emit:
        reports = filter_info_by_target_chr2.out
}

// ================================================================================
// MAIN WORKFLOW
// ================================================================================

workflow {
    // Display pipeline header
    displayHeader()
    
    // ────────────────────────────────────────
    // STAGE 1: Data Preparation
    // ────────────────────────────────────────
    preprocess(params.target_datasets)
    
    // ────────────────────────────────────────
    // STAGE 2: Data Chunking
    // ────────────────────────────────────────
    subset(preprocess.out.dataset_qc)
    
    // ────────────────────────────────────────
    // STAGE 3: Phasing
    // ────────────────────────────────────────
    phasing_data = Channel.from(params.ref_panels)
        .combine(subset.out.chunks)
        .flatMap{ ref_name, ref_m3vcf, ref_vcf, dataset, chrm, start, end, dataset_, dataset_vcf ->
            def vcf = sprintf(ref_vcf, chrm)
            def m3vcf = sprintf(ref_m3vcf, chrm)
            def vcf_idx = vcf.endsWith("vcf.gz") ? "${vcf}.tbi" : "${vcf}.csi"
            
            [[chrm, ref_name, file(m3vcf), file(vcf), file(vcf_idx), 
              file(params.eagle_genetic_map), start, end, dataset, dataset, file(dataset_vcf)]]
        }
    
    phasing(phasing_data)
    
    // ────────────────────────────────────────
    // STAGE 4: Imputation
    // ────────────────────────────────────────
    impute(phasing.out.chunks_phased)
    
    // ────────────────────────────────────────
    // STAGE 5: Reporting
    // ────────────────────────────────────────
    
    // Prepare imputed data for reporting
    impute_data = impute.out.chunks_imputed
        .map{ chr, fwd, rev, test_data, ref, imputed_vcf, imputed_info, tst_data -> 
            [test_data, ref, imputed_vcf, imputed_info]
        }
        .combine(params.target_datasets, by: 0)
        .map{ test_data, ref, imputed_vcf, imputed_info, orig_vcf -> 
            [test_data, ref, orig_vcf, imputed_vcf, imputed_info]
        }
    
    // Genome-wide reports
    report_by_ref(impute_data)
    report_by_dataset(impute_data)
    
    // Chromosome-level reports
    impute_data_chr = impute.out.chunks_imputed
        .map{ chrm, chunk_start, chunk_end, target_name, ref_name, imputed_bcf, info_file, tagName -> 
            [target_name, chrm, ref_name, imputed_bcf, info_file]
        }
        .combine(Channel.from(params.target_datasets), by: 0)
        .map{ target_name, chrm, ref_name, imputed_vcf, imputed_info, orig_vcf -> 
            [chrm, target_name, ref_name, orig_vcf, imputed_vcf, imputed_info]
        }
    
    report_by_ref_chromosome(impute_data_chr)
    report_by_dataset_chromosome(impute_data_chr)
    
    // ────────────────────────────────────────
    // STAGE 6: Frequency Analysis
    // ────────────────────────────────────────
    impute_data_with_chr = impute.out.chunks_imputed
        .map{ chr, fwd, rev, test_data, ref, imputed_vcf, imputed_info, tst_data -> 
            [chr, test_data, ref, imputed_vcf, imputed_info]
        }
    
    // Generate frequency data
    freq_input = impute_data_with_chr
        .map{ chr, target_name, ref_name, impute_vcf, info -> 
            def ref_panel = params.ref_panels.find { it[0] == ref_name }
            if (ref_panel) {
                def ref_vcf_path = ref_panel[2].replace('%s', chr)
                [target_name, ref_name, file(impute_vcf), file(ref_vcf_path)]
            }
        }
        .filter { it != null }
    
    generate_frequency(freq_input)
    
    // Frequency comparison plots
    freq_comp = impute_data_with_chr
        .map{ chr, target_name, ref_name, impute_vcf, info -> 
            [target_name, ref_name, info]
        }
        .combine(generate_frequency.out, by: [0,1])
    
    plot_freq_comparison(freq_comp)
    
    // R2 vs SNP position plots
    combineInfo_frq = impute_data_with_chr
        .map{ chr, target_name, ref_name, impute_vcf, info ->
            [target_name, ref_name, info, params.maf_thresh]
        }
        .combine(generate_frequency.out, by: [0,1])
        .map{ target_name, ref_name, info, maf_thresh, target_frq, ref_frq -> 
            [target_name, ref_name, info, maf_thresh, target_frq]
        }
    
    plot_r2_SNPpos(combineInfo_frq)
    
    // Average R-squared calculation
    rsquared_input = impute_data
        .map{ target_name, ref_name, vcf, impute_vcf, info ->
            [target_name, ref_name, info]
        }
    
    average_r2(rsquared_input)
    
    log.info ""
    log.info "═════════════════════════════════════════════════════════════"
    log.info "Pipeline execution completed successfully!"
    log.info "Results available in: ${params.outDir}"
    log.info "═════════════════════════════════════════════════════════════"
}