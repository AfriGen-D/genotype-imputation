/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORT_DATASET
    Dataset-level aggregation of chunk reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COMBINE_CHUNK_REPORTS     } from '../../modules/local/report/combine_chunk_reports'
include { SUMMARIZE_DATASET_METRICS } from '../../modules/local/report/summarize_dataset_metrics'
include { PLOT_DATASET_SUMMARY      } from '../../modules/local/report/plot_dataset_summary'

workflow REPORT_DATASET {
    take:
    ch_chunk_reports // channel: [ val(meta), val(ref_name), path(reports) ]
    ch_chunk_plots   // channel: [ val(meta), val(ref_name), path(plots) ]

    main:
    ch_versions = Channel.empty()
    
    //
    // Group chunk reports by dataset (using meta.sample as the dataset identifier)
    //
    ch_reports_by_dataset = ch_chunk_reports
        .map { meta, ref_name, reports ->
            // Create new meta with dataset-level ID
            def dataset_meta = [:]
            dataset_meta.id = meta.sample ?: meta.id.split('_')[0]  // Extract dataset name
            dataset_meta.sample = meta.sample ?: meta.id.split('_')[0]
            dataset_meta.population = meta.population
            dataset_meta.study = meta.study
            [dataset_meta, ref_name, reports]
        }
        .groupTuple(by: [0, 1])  // Group by dataset_meta and ref_name
    
    ch_plots_by_dataset = ch_chunk_plots
        .map { meta, ref_name, plots ->
            // Create new meta with dataset-level ID
            def dataset_meta = [:]
            dataset_meta.id = meta.sample ?: meta.id.split('_')[0]  // Extract dataset name
            dataset_meta.sample = meta.sample ?: meta.id.split('_')[0]
            dataset_meta.population = meta.population
            dataset_meta.study = meta.study
            [dataset_meta, ref_name, plots]
        }
        .groupTuple(by: [0, 1])  // Group by dataset_meta and ref_name
    
    //
    // MODULE: Combine chunk reports into dataset-level report
    //
    COMBINE_CHUNK_REPORTS ( ch_reports_by_dataset )
    ch_versions = ch_versions.mix(COMBINE_CHUNK_REPORTS.out.versions)
    
    //
    // MODULE: Generate dataset-level summary metrics
    //
    SUMMARIZE_DATASET_METRICS ( COMBINE_CHUNK_REPORTS.out.combined )
    ch_versions = ch_versions.mix(SUMMARIZE_DATASET_METRICS.out.versions)
    
    //
    // MODULE: Create dataset-level summary plots
    //
    ch_summary_input = SUMMARIZE_DATASET_METRICS.out.summary
        .join(ch_plots_by_dataset, by: [0, 1], remainder: true)
        .map { meta, ref_name, summary, plots ->
            [meta, ref_name, summary, plots ?: []]
        }
    
    PLOT_DATASET_SUMMARY ( ch_summary_input )
    ch_versions = ch_versions.mix(PLOT_DATASET_SUMMARY.out.versions)
    
    emit:
    dataset_reports  = COMBINE_CHUNK_REPORTS.out.combined
    dataset_summary  = SUMMARIZE_DATASET_METRICS.out.summary
    dataset_plots    = PLOT_DATASET_SUMMARY.out.plot
    versions        = ch_versions
}