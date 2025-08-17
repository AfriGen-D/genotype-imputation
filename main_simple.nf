#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    h3abionet/chipimputation - Simplified pipeline with inline workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import QC modules
include { get_chromosome; fill_tags_vcf; check_chromosome; check_files; check_chromosome_vcf; 
         check_mismatch; no_mismatch; qc_dupl; split_multi_allelic; filter_min_ac; 
         target_qc; qc_site_missingness; sites_only; combine_vcfs; combine_infos; 
         combine_csvs as combine_freqs; combine_vcfs_chrm } from './modules/qc'

// Import VCF processing modules         
include { vcf_map_simple; extract_site_from_vcf; generate_chunks_vcf; split_target_to_chunk; 
         vcf_map; vcf_freq; info_freq; fill_tags_VCF; sort_vcf; get_vcf_sites; 
         extract_pop } from './modules/subset_vcf'

// Import phasing modules         
include { minimac4_phasing_eagle } from './modules/phasing'

// Import imputation modules
include { impute_minimac4; extract_impute_info; impute_minimac4_1; combineImpute; 
         combineInfo; filter_info_by_target } from './modules/impute'

// Import reporting modules         
include { filter_info; report_site_by_maf; plot_freq_comparison; report_well_imputed_by_target; 
         plot_performance_target; report_accuracy_target; plot_accuracy_target; 
         generate_frequency; plot_r2_SNPpos; plot_r2_SNPcount; plot_hist_r2_SNPcount; 
         plot_MAF_r2; average_r2 } from './modules/report'
         
include { filter_info_by_target_chr; filter_info_by_target_chr2; report_well_imputed_by_target_chr; 
         report_well_imputed_by_target_chr2; plot_performance_target_chr; plot_performance_target_chr2; 
         report_accuracy_target_chr; report_accuracy_target_chr2; plot_accuracy_target_chr; 
         plot_accuracy_target_chr2; plot_r2_SNPcount_chr; plot_hist_r2_SNPcount_chr; 
         plot_MAF_r2_chr } from './modules/report_chr'
         
include { run_qc_plots } from './modules/qc_plots'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER DEFAULTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Set parameter defaults
params.help                 = false
params.version             = '1.0.0'
params.input               = null
params.outDir              = './results'
params.project_name        = 'imputation_project'

// Data processing parameters
params.chromosomes         = 'ALL'
params.chunk_size          = 20000000  // 20MB chunks
params.buffer_size         = 1000000   // 1MB buffer

// QC parameters
params.site_miss           = 0.05
params.hwe                 = 0.00001
params.mac                 = 1
params.min_ac              = 2
params.min_alleles         = 2
params.max_alleles         = 2
params.maf_thresh          = 0.1

// Phasing parameters
params.phasing_method      = 'eagle'
params.eagle_pbwt_iters    = 2
params.eagle_genetic_map   = null

// Imputation parameters
params.impute_method       = 'minimac4'
params.minRatio           = 0.01
params.NE                 = 20000
params.impute_iter        = 10
params.impute_burnin      = 2
params.impute_info_cutoff = 0.3

// Reference data
params.reference_genome    = null
params.ref_panels         = []
params.target_datasets    = []

// Reporting parameters
params.qc_plots           = true
params.r2_threshold       = 0.3

// Compute resources
params.max_cpus           = 16
params.max_memory         = '128.GB'
params.max_time           = '240.h'

// Other
params.email              = null
params.plink              = 'plink2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info """
    =========================================
    h3abionet/chipimputation v${params.version}
    =========================================
    
    Usage:
    nextflow run main_simple.nf --input samplesheet.csv --outDir results -profile docker
    
    OR with config file:
    nextflow run main_simple.nf -c myconfig.config -profile slurm,singularity
    
    Mandatory arguments:
      --input               Path to input samplesheet CSV file OR use --target_datasets
      --outDir              The output directory where results will be saved
      --eagle_genetic_map   Path to Eagle genetic map file
      --reference_genome    Path to reference genome FASTA file
      --ref_panels          Reference panels for imputation (list format)
      
    See full help with --help
    """.stripIndent()
    System.exit(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Function to check if files exist
def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}

// Validate input
if (!params.input && !params.target_datasets) {
    log.error "ERROR: Please provide either --input (samplesheet) or --target_datasets"
    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW DEFINITIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow preprocess {
    take: datasets

    main:
        // Check if study genotype files exist
        target_datasets = []
        datasets.each { name, vcf ->
            check_files([vcf])
            target_datasets << [name, file(vcf)]
        }
        target_datasets = Channel.from(target_datasets)

        // Check if eagle map file exists
        if(params.eagle_genetic_map) {
            check_files([params.eagle_genetic_map])
        }

        // Check if fasta reference genome exists
        if(params.reference_genome) {
            check_files([params.reference_genome, "${params.reference_genome}.fai"])
        }

        // Check reference panel files
        def resolved_chromosomes = []
        if (params.chromosomes == 'ALL' || params.chromosomes == '') {
            log.info "|-- INFO: Chromosomes set to 'ALL' - checking reference files"
            
            def b37_chrms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
            def b38_chrms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
            
            def test_chrm_found = false
            params.ref_panels.each{ ref_name, ref_msav, ref_vcf ->
                if (!test_chrm_found) {
                    def test_vcf_b38 = sprintf(ref_vcf, 'chr1')
                    def test_msav_b38 = sprintf(ref_msav, 'chr1')
                    if (file(test_vcf_b38).exists() || file(test_msav_b38).exists()) {
                        resolved_chromosomes = b38_chrms
                        test_chrm_found = true
                        log.info "|-- INFO: Detected b38 chromosome format"
                    } else {
                        def test_vcf_b37 = sprintf(ref_vcf, '1')
                        def test_msav_b37 = sprintf(ref_msav, '1')
                        if (file(test_vcf_b37).exists() || file(test_msav_b37).exists()) {
                            resolved_chromosomes = b37_chrms
                            test_chrm_found = true
                            log.info "|-- INFO: Detected b37 chromosome format"
                        }
                    }
                }
            }
            
            if (!test_chrm_found) {
                resolved_chromosomes = b37_chrms + b38_chrms
                log.warn "|-- WARN: Could not determine chromosome format"
            }
        } else {
            resolved_chromosomes = params.chromosomes.split(',')
        }

        // QC steps
        check_mismatch(target_datasets.map{ dataset, dataset_vcf -> 
            [ dataset, '', '', '', file(dataset_vcf), file(params.reference_genome) ] 
        })
        
        qc_dupl(target_datasets.map{ dataset, dataset_vcf -> 
            [ dataset, '', '', '', file(dataset_vcf) ] 
        })
        
        split_multi_allelic(qc_dupl.out)
        fill_tags_vcf(split_multi_allelic.out)
        filter_min_ac(fill_tags_vcf.out.map{ dataset, chrm, start, end, vcf -> 
            [ dataset, chrm, start, end, file(vcf), 
              " --min-ac ${params.min_ac} --max-alleles ${params.max_alleles} --min-alleles ${params.min_alleles} -v snps " ] 
        })

    emit:
        dataset_qc = filter_min_ac.out
}

workflow subset {
    take: data
    
    main:
        get_chromosome(data.map{ dataset, chrm, start, end, vcf, map_file -> 
            [ dataset, file(vcf) ] 
        })
        
        data_chrms = get_chromosome.out
            .map{ dataset, dataset_vcf, map_file -> 
                check_chromosome_vcf(dataset, dataset_vcf, map_file, params.chromosomes) 
            }
            .map{ dataset, dataset_vcf, map_file, chrms -> 
                [ dataset, file(dataset_vcf), file(map_file), chrms.unique().join(',') ] 
            }
        
        generate_chunks_vcf(data_chrms.map{ dataset, vcf, map_file, chrms -> 
            [ dataset, file(vcf), file(map_file), chrms, params.chunk_size ] 
        })
        
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

workflow phasing {
    take: data
    
    main:
        data.map{ dataset, chrm, start, end, dataset_vcf ->
            ref_panels = []
            params.ref_panels.each{ ref_name, ref_msav, ref_vcf ->
                ref_panels << [dataset, chrm, start, end, ref_name, file(dataset_vcf), 
                              file(sprintf(ref_vcf, chrm)), file(sprintf(ref_msav, chrm))]
            }
            return ref_panels
        }.flatMap().set{ data_ref }
        
        minimac4_phasing_eagle(data_ref)
    
    emit: 
        phased_data = minimac4_phasing_eagle.out
}

workflow impute {
    take: phased_data 
    
    main:
        impute_minimac4(phased_data)
        extract_impute_info(impute_minimac4.out)
    
    emit:
        imputed_data = impute_minimac4.out
        impute_info = extract_impute_info.out
}

workflow report_by_ref {
    take:
        imputed_data
        impute_info
        dataset_qc
    
    main:
        // Group by reference panel
        imputeCombine_ref = imputed_data
            .map{ dataset, chrm, start, end, ref_name, imputed_vcf, info ->
                [ref_name, dataset, imputed_vcf, info]
            }
            .groupTuple(by: 0)
            .map{ ref_name, datasets, vcfs, infos ->
                [ref_name, datasets.join(','), '', infos.join(',')]
            }
        
        filter_info_by_target(imputeCombine_ref)
        
        report_well_imputed_by_target(filter_info_by_target.out.map{ 
            target_name, ref_panels, wellInfo, accInfo -> 
            [ target_name, ref_panels, file(wellInfo) ]
        })
        
        plot_performance_target(report_well_imputed_by_target.out.map{ 
            target_name, ref_panels, wellInfo, wellInfo_summary -> 
            [ target_name, ref_panels, file(wellInfo), file(wellInfo_summary), 'DATASETS' ]
        })
        
        report_accuracy_target(filter_info_by_target.out.map{ 
            target_name, ref_panels, wellInfo, accInfo -> 
            [ target_name, ref_panels, file(accInfo), 'DATASETS' ]
        })
        
        plot_accuracy_target(report_accuracy_target.out)
        
    emit:
        reports = filter_info_by_target.out
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def intro() {
    log.info ""
    log.info """
    =======================================================
    h3achipimputation v${params.version}
    ======================================================= """
    def summary = [:]
    summary['Pipeline Name']    = 'h3achipimputation'
    summary['Pipeline version'] = params.version
    summary['Run Name']         = workflow.runName
    
    if(params.target_datasets) {
        summary['Target datasets'] = params.target_datasets.collect{ "${it[0]} (${it[1]})" }.join(', ')
    } else if(params.input) {
        summary['Input'] = params.input
    }
    
    if(params.ref_panels) {
        summary['Reference panels'] = params.ref_panels.collect{ "${it[0]}" }.join(', ')
    }
    
    summary['Output dir']       = params.outDir
    summary['Chromosomes']      = params.chromosomes
    summary['Chunk Size']       = params.chunk_size
    summary['Max Memory']       = params.max_memory
    summary['Max CPUs']         = params.max_cpus
    summary['Max Time']         = params.max_time
    summary['Working dir']      = workflow.workDir
    summary['Container Engine'] = workflow.containerEngine ?: 'none'
    summary['Config Profile']   = workflow.profile
    
    log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
    log.info "======================================================="
    log.info ""
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    intro()
    
    // Create input channel
    if (params.input) {
        // Read from CSV samplesheet
        Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                [ row.sample, file(row.vcf) ]
            }
            .collect()
            .set { ch_input_list }
    } else {
        // Use target_datasets parameter directly
        ch_input_list = params.target_datasets
    }
    
    // Step 1: Preprocessing and QC
    preprocess(ch_input_list)
    
    // Step 2: Subset data into chunks
    subset(preprocess.out.dataset_qc)
    
    // Step 3: Phasing
    phasing(subset.out.chunks)
    
    // Step 4: Imputation
    impute(phasing.out.phased_data)
    
    // Step 5: Reporting and QC plots
    if (params.qc_plots) {
        report_by_ref(
            impute.out.imputed_data,
            impute.out.impute_info,
            preprocess.out.dataset_qc
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION NOTIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    log.info ""
    log.info "======================================================="
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Duration    : ${workflow.duration}"
    log.info "Success     : ${workflow.success}"
    log.info "Exit status : ${workflow.exitStatus}"
    log.info "Output dir  : ${params.outDir}"
    log.info ""
    
    if (workflow.success) {
        log.info "Pipeline completed successfully!"
    } else {
        log.info "Pipeline completed with errors"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/