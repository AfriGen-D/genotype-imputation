#!/usr/bin/env nextflow
/*
 * Test script to verify PLOT_FREQ_COMPARISON module works
 */

nextflow.enable.dsl = 2

// Import the module
include { PLOT_FREQ_COMPARISON } from './modules/local/report/plot_freq_comparison'

// Test workflow
workflow {
    // Create test channel with frequency comparison TSV file
    ch_test = Channel.of(
        [
            [id: 'test_quick_chr21'],
            'chunk',
            file('/scratch3/users/mamana/nextflow-work/c5/40cd81596249307cacc788bc1767eb/test_quick_chr21_15000001_15030080_H3AR6x_chr21_15000001_15030080.freq_comparison.tsv')
        ]
    )
    
    // Run the module
    PLOT_FREQ_COMPARISON(ch_test)
    
    // View output
    PLOT_FREQ_COMPARISON.out.plot.view()
}