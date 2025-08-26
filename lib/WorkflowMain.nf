/*
========================================================================================
    WORKFLOW HELPER FUNCTIONS
========================================================================================
*/

class WorkflowMain {

    //
    // Check and validate parameters
    //
    public static void initialise(workflow, params, log) {
        
        if (!params.input) {
            log.error "Please provide an input samplesheet with the --input parameter"
            System.exit(1)
        }

        if (!params.reference_panels || params.reference_panels.size() == 0) {
            log.error "Please provide reference panels with the --reference_panels parameter"
            System.exit(1)
        }

        // Validate genome build
        if (params.genome_build && !(params.genome_build in ['b37', 'b38', 'hg19', 'hg38'])) {
            log.error "ERROR: Invalid genome build '${params.genome_build}'. Must be 'b37', 'b38', 'hg19', or 'hg38'"
            System.exit(1)
        }
        
        // Validate imputation tool
        if (params.imputation_tool && !(params.imputation_tool in ['minimac4', 'impute5', 'beagle5'])) {
            log.error "ERROR: Invalid imputation tool '${params.imputation_tool}'. Must be 'minimac4', 'impute5', or 'beagle5'"
            System.exit(1)
        }
        
        // Validate phasing tool
        if (params.phasing_tool && !(params.phasing_tool in ['eagle', 'shapeit', 'beagle'])) {
            log.error "ERROR: Invalid phasing tool '${params.phasing_tool}'. Must be 'eagle', 'shapeit', or 'beagle'"
            System.exit(1)
        }

        // Print workflow summary
        log.info header(workflow)
        log.info summary(workflow, params)
    }

    //
    // Get workflow summary header
    //
    public static String header(workflow) {
        def header_text = """
        ====================================
         ${workflow.manifest.name} v${workflow.manifest.version}
        ====================================
        """.stripIndent()
        return header_text
    }

    //
    // Get workflow summary
    //
    public static String summary(workflow, params) {
        def summary_text = """
        Input            : ${params.input}
        Output directory : ${params.outdir}
        Genome build     : ${params.genome_build}
        Phasing tool     : ${params.phasing_tool}
        Imputation tool  : ${params.imputation_tool}
        Report level     : ${params.report_level ?: 'standard'}
        Max CPUs         : ${params.max_cpus}
        Max Memory       : ${params.max_memory}
        Container        : ${workflow.containerEngine ?: 'None'}
        Profile          : ${workflow.profile}
        """.stripIndent()
        return summary_text
    }

    //
    // Get attribute from genome config
    //
    public static String getGenomeAttribute(params, attribute) {
        if (!params.genomes) {
            return params[attribute]
        }
        if (!params.genome) {
            return params[attribute]
        }
        if (!params.genomes.containsKey(params.genome)) {
            return params[attribute]
        }
        return params.genomes[params.genome][attribute]
    }
}