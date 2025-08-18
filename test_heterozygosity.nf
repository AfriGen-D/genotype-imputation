#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PLOT_HETEROZYGOSITY } from './modules/local/report/plot_heterozygosity'

workflow {
    println "Testing plot_heterozygosity module"
}