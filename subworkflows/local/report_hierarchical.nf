/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORT_AGG
    Hierarchical aggregation: Chunk → Chromosome → Genome-wide
    All aggregation happens AFTER imputation is complete
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AGGREGATE_CHUNKS_TO_CHR   } from '../../modules/local/report/aggregate_chunks_to_chr'
include { AGGREGATE_CHR_TO_GENOME   } from '../../modules/local/report/aggregate_chr_to_genome'
include { PLOT_CHR_SUMMARY          } from '../../modules/local/report/plot_chr_summary'
include { PLOT_CHR_ALL_METRICS      } from '../../modules/local/report/plot_chr_all_metrics'
include { PLOT_GENOME_SUMMARY       } from '../../modules/local/report/plot_genome_summary'
include { PLOT_GENOME_ALL_METRICS   } from '../../modules/local/report/plot_genome_all_metrics'
include { GENERATE_DATASET_REPORT   } from '../../modules/local/report/generate_dataset_report'
include { GENERATE_REPORT_INDEX     } from '../../modules/local/report/generate_report_index'

workflow REPORT_AGG {
    take:
    ch_chunk_reports // channel: [ val(meta), val(ref_name), path(reports) ] - post-imputation reports
    ch_chunk_plots   // channel: [ val(meta), val(ref_name), path(plots) ]   - post-imputation plots
    ch_chunk_info    // channel: [ val(meta), val(ref_name), path(info) ]    - imputation info files

    main:
    ch_versions = Channel.empty()
    
    //
    // STEP 1: Group chunks by chromosome within each dataset
    // meta.id format: dataset_chr##_start_end
    //
    ch_reports_by_chr = ch_chunk_reports
        .map { meta, ref_name, reports ->
            def dataset = meta.sample ?: meta.id.split('_chr')[0]
            def chr_match = meta.id =~ /chr(\d+)/
            def chromosome = chr_match ? "chr${chr_match[0][1]}" : "unknown"
            
            def chr_meta = [:]
            chr_meta.id = "${dataset}_${chromosome}"
            chr_meta.dataset = dataset
            chr_meta.chromosome = chromosome
            chr_meta.population = meta.population
            chr_meta.study = meta.study
            
            [chr_meta, ref_name, reports]
        }
        .groupTuple(by: [0, 1])  // Group by chr_meta and ref_name
        .map { chr_meta, ref_name, report_lists ->
            [chr_meta, ref_name, report_lists.flatten()]
        }
    
    ch_plots_by_chr = ch_chunk_plots
        .map { meta, ref_name, plots ->
            def dataset = meta.sample ?: meta.id.split('_chr')[0]
            def chr_match = meta.id =~ /chr(\d+)/
            def chromosome = chr_match ? "chr${chr_match[0][1]}" : "unknown"
            
            def chr_meta = [:]
            chr_meta.id = "${dataset}_${chromosome}"
            chr_meta.dataset = dataset
            chr_meta.chromosome = chromosome
            chr_meta.population = meta.population
            chr_meta.study = meta.study
            
            [chr_meta, ref_name, plots]
        }
        .groupTuple(by: [0, 1])
        .map { chr_meta, ref_name, plot_lists ->
            [chr_meta, ref_name, plot_lists.flatten()]
        }
    
    ch_info_by_chr = ch_chunk_info
        .map { meta, ref_name, info ->
            def dataset = meta.sample ?: meta.id.split('_chr')[0]
            def chr_match = meta.id =~ /chr(\d+)/
            def chromosome = chr_match ? "chr${chr_match[0][1]}" : "unknown"
            
            def chr_meta = [:]
            chr_meta.id = "${dataset}_${chromosome}"
            chr_meta.dataset = dataset
            chr_meta.chromosome = chromosome
            chr_meta.population = meta.population
            chr_meta.study = meta.study
            
            [chr_meta, ref_name, info]
        }
        .groupTuple(by: [0, 1])
        .map { chr_meta, ref_name, info_lists ->
            [chr_meta, ref_name, info_lists.flatten()]
        }
    
    //
    // MODULE: Aggregate chunks to chromosome level
    //
    AGGREGATE_CHUNKS_TO_CHR ( 
        ch_reports_by_chr.join(ch_info_by_chr, by: [0, 1])
    )
    ch_versions = ch_versions.mix(AGGREGATE_CHUNKS_TO_CHR.out.versions)
    
    //
    // MODULE: Create chromosome-level plots
    //
    // Use leftJoin to handle empty chunk plots when they're disabled
    ch_chr_summary_with_plots = AGGREGATE_CHUNKS_TO_CHR.out.chr_summary
        .join(ch_plots_by_chr, by: [0, 1], remainder: true)
        .map { chr_meta, ref_name, chr_summary, plots ->
            // If plots is null (chunk plots disabled), provide empty list
            [chr_meta, ref_name, chr_summary, plots ?: []]
        }
    
    PLOT_CHR_SUMMARY ( ch_chr_summary_with_plots )
    ch_versions = ch_versions.mix(PLOT_CHR_SUMMARY.out.versions)
    
    //
    // MODULE: Create comprehensive chromosome-level plots (all metrics)
    //
    PLOT_CHR_ALL_METRICS ( ch_chr_summary_with_plots )
    ch_versions = ch_versions.mix(PLOT_CHR_ALL_METRICS.out.versions)
    
    //
    // STEP 2: Group chromosomes by dataset for genome-wide analysis
    //
    ch_chr_reports_by_dataset = AGGREGATE_CHUNKS_TO_CHR.out.chr_summary
        .map { chr_meta, ref_name, chr_summary ->
            def genome_meta = [:]
            genome_meta.id = chr_meta.dataset
            genome_meta.dataset = chr_meta.dataset
            genome_meta.population = chr_meta.population
            genome_meta.study = chr_meta.study
            
            [genome_meta, ref_name, chr_meta.chromosome, chr_summary]
        }
        .groupTuple(by: [0, 1])  // Group by genome_meta and ref_name
        .map { genome_meta, ref_name, chromosomes, summaries ->
            [genome_meta, ref_name, summaries.flatten()]
        }
    
    ch_chr_plots_by_dataset = PLOT_CHR_SUMMARY.out.chr_plots
        .map { chr_meta, ref_name, chr_plots ->
            def genome_meta = [:]
            genome_meta.id = chr_meta.dataset
            genome_meta.dataset = chr_meta.dataset
            genome_meta.population = chr_meta.population
            genome_meta.study = chr_meta.study
            
            [genome_meta, ref_name, chr_plots]
        }
        .groupTuple(by: [0, 1])
        .map { genome_meta, ref_name, plot_lists ->
            [genome_meta, ref_name, plot_lists.flatten()]
        }
    
    //
    // MODULE: Aggregate chromosome data to genome-wide level
    //
    AGGREGATE_CHR_TO_GENOME ( ch_chr_reports_by_dataset )
    ch_versions = ch_versions.mix(AGGREGATE_CHR_TO_GENOME.out.versions)
    
    //
    // MODULE: Create genome-wide plots
    //
    // Use leftJoin to handle empty chromosome plots
    ch_genome_summary_with_plots = AGGREGATE_CHR_TO_GENOME.out.genome_summary
        .join(ch_chr_plots_by_dataset, by: [0, 1], remainder: true)
        .map { genome_meta, ref_name, genome_summary, plots ->
            // If plots is null, provide empty list
            [genome_meta, ref_name, genome_summary, plots ?: []]
        }
    
    PLOT_GENOME_SUMMARY ( ch_genome_summary_with_plots )
    ch_versions = ch_versions.mix(PLOT_GENOME_SUMMARY.out.versions)
    
    //
    // MODULE: Create comprehensive genome-wide plots (all metrics)
    //
    PLOT_GENOME_ALL_METRICS ( ch_genome_summary_with_plots )
    ch_versions = ch_versions.mix(PLOT_GENOME_ALL_METRICS.out.versions)
    
    //
    // MODULE: Generate comprehensive dataset report (PDF/HTML)
    //
    ch_final_report_input = AGGREGATE_CHR_TO_GENOME.out.genome_summary
        .join(PLOT_GENOME_SUMMARY.out.genome_plots, by: [0, 1])
    
    GENERATE_DATASET_REPORT ( ch_final_report_input )
    ch_versions = ch_versions.mix(GENERATE_DATASET_REPORT.out.versions)
    
    //
    // Generate HTML index for all reports
    //
    ch_dataset_for_index = ch_final_report_input
        .map { meta, ref_name, genome_summary, genome_plots -> 
            meta.dataset ?: meta.id
        }
        .unique()
        .take(1)
    
    ch_reports_dir = Channel.value(file("${params.outdir}/reports"))
    
    GENERATE_REPORT_INDEX (
        ch_dataset_for_index,
        ch_reports_dir
    )
    ch_versions = ch_versions.mix(GENERATE_REPORT_INDEX.out.versions)
    
    emit:
    chr_summaries        = AGGREGATE_CHUNKS_TO_CHR.out.chr_summary
    chr_plots           = PLOT_CHR_SUMMARY.out.chr_plots
    chr_all_plots       = PLOT_CHR_ALL_METRICS.out.combined
    chr_performance     = PLOT_CHR_ALL_METRICS.out.performance
    chr_r2_position     = PLOT_CHR_ALL_METRICS.out.r2_position
    chr_accuracy        = PLOT_CHR_ALL_METRICS.out.accuracy
    chr_maf_analysis    = PLOT_CHR_ALL_METRICS.out.maf_analysis
    chr_freq_comparison = PLOT_CHR_ALL_METRICS.out.freq_comparison
    genome_summary      = AGGREGATE_CHR_TO_GENOME.out.genome_summary
    genome_plots        = PLOT_GENOME_SUMMARY.out.genome_plots
    genome_all_plots    = PLOT_GENOME_ALL_METRICS.out.combined
    genome_performance  = PLOT_GENOME_ALL_METRICS.out.performance
    genome_r2_dist      = PLOT_GENOME_ALL_METRICS.out.r2_distribution
    genome_accuracy     = PLOT_GENOME_ALL_METRICS.out.accuracy
    genome_maf_analysis = PLOT_GENOME_ALL_METRICS.out.maf_analysis
    genome_freq_comp    = PLOT_GENOME_ALL_METRICS.out.freq_comparison
    genome_chr_comp     = PLOT_GENOME_ALL_METRICS.out.chr_comparison
    final_report        = GENERATE_DATASET_REPORT.out.report
    report_index        = GENERATE_REPORT_INDEX.out.main_index
    versions            = ch_versions
}