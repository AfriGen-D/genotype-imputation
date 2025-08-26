/*
========================================================================================
    MAIN WORKFLOW HELPER FUNCTIONS
========================================================================================
*/

class WorkflowMain {

    /*
     * Validate parameters and print summary
     */
    public static String validateParameters(workflow, params, log) {
        
        // Check mandatory parameters
        if (!params.input) {
            log.error "ERROR: Input samplesheet not specified with --input"
            System.exit(1)
        }
        
        if (!params.reference_panels || params.reference_panels.size() == 0) {
            log.error "ERROR: Reference panels not specified with --reference_panels"
            System.exit(1)
        }
        
        // Validate genome build
        if (!(params.genome_build in ['b37', 'b38'])) {
            log.error "ERROR: Invalid genome build '${params.genome_build}'. Must be 'b37' or 'b38'"
            System.exit(1)
        }
        
        // Validate tools
        def valid_phasing_tools = ['eagle', 'shapeit', 'beagle']
        if (!(params.phasing_tool in valid_phasing_tools)) {
            log.error "ERROR: Invalid phasing tool '${params.phasing_tool}'. Must be one of: ${valid_phasing_tools.join(', ')}"
            System.exit(1)
        }
        
        def valid_imputation_tools = ['minimac4', 'impute5', 'beagle5']
        if (!(params.imputation_tool in valid_imputation_tools)) {
            log.error "ERROR: Invalid imputation tool '${params.imputation_tool}'. Must be one of: ${valid_imputation_tools.join(', ')}"
            System.exit(1)
        }
        
        // Validate report level
        def valid_report_levels = ['summary', 'detailed', 'full']
        if (!(params.report_level in valid_report_levels)) {
            log.error "ERROR: Invalid report level '${params.report_level}'. Must be one of: ${valid_report_levels.join(', ')}"
            System.exit(1)
        }
        
        return "Parameters validated successfully"
    }

    /*
     * Generate parameter summary for logging
     */
    public static String paramsSummaryLog(workflow, params, log) {
        def summary = [:]
        summary['Pipeline Name']    = workflow.manifest.name
        summary['Pipeline Version'] = workflow.manifest.version
        summary['Run Name']         = workflow.runName
        summary['Input']            = params.input
        summary['Genome Build']     = params.genome_build
        summary['Reference Panels'] = params.reference_panels.collect{ it[0] }.join(', ')
        summary['Phasing Tool']     = params.phasing_tool
        summary['Imputation Tool']  = params.imputation_tool
        summary['Max CPUs']         = params.max_cpus
        summary['Max Memory']       = params.max_memory
        summary['Max Time']         = params.max_time
        summary['Output dir']       = params.outdir
        summary['Work dir']         = workflow.workDir
        summary['Container Engine'] = workflow.containerEngine ?: 'None'
        summary['Config Profile']   = workflow.profile
        summary['User']             = workflow.userName
        
        def summary_log = []
        summary_log << "\\033[0;35m========================================\\033[0m"
        summary_log << "\\033[0;35m${workflow.manifest.name} v${workflow.manifest.version}\\033[0m"
        summary_log << "\\033[0;35m========================================\\033[0m"
        summary.each{ k, v ->
            def key = k.padRight(20)
            summary_log << "${key}: ${v}"
        }
        summary_log << "\\033[0;35m========================================\\033[0m"
        
        return summary_log.join("\\n")
    }

    /*
     * Parse samplesheet and return channel
     */
    public static LinkedHashMap parseSamplesheet(samplesheet_file) {
        def samplesheet = []
        def header = []
        def line_count = 0
        
        samplesheet_file.eachLine { line ->
            line_count++
            def row = line.tokenize(',')
            
            if (line_count == 1) {
                header = row
                // Validate required columns
                def required = ['sample', 'vcf']
                def missing = required - header
                if (missing.size() > 0) {
                    throw new Exception("ERROR: Missing required columns in samplesheet: ${missing.join(', ')}")
                }
            } else {
                def meta = [:]
                row.eachWithIndex { val, idx ->
                    meta[header[idx]] = val
                }
                
                // Validate file exists
                if (!file(meta.vcf).exists()) {
                    throw new Exception("ERROR: VCF file does not exist: ${meta.vcf}")
                }
                
                // Add metadata
                meta.id = meta.sample
                meta.population = meta.population ?: 'Unknown'
                meta.sex = meta.sex ?: 'Unknown'
                
                samplesheet << meta
            }
        }
        
        return samplesheet
    }

    /*
     * Check file exists
     */
    public static void checkFileExists(file_path, param_name, log) {
        if (!file(file_path).exists()) {
            log.error "ERROR: File specified with ${param_name} does not exist: ${file_path}"
            System.exit(1)
        }
    }

    /*
     * Generate completion email HTML
     */
    public static String completionEmail(workflow, params, summary_log, log) {
        
        def status = workflow.success ? 'SUCCESS' : 'FAILED'
        def color = workflow.success ? 'green' : 'red'
        
        def html = """
        <html>
        <head>
            <style>
                body { font-family: Arial, sans-serif; }
                h2 { color: ${color}; }
                table { border-collapse: collapse; width: 100%; }
                th, td { text-align: left; padding: 8px; border-bottom: 1px solid #ddd; }
                th { background-color: #f2f2f2; }
                .success { color: green; }
                .failed { color: red; }
            </style>
        </head>
        <body>
            <h2>Imputation Pipeline: ${status}</h2>
            <h3>Run Summary</h3>
            <table>
                <tr><th>Parameter</th><th>Value</th></tr>
                <tr><td>Pipeline</td><td>${workflow.manifest.name} v${workflow.manifest.version}</td></tr>
                <tr><td>Run Name</td><td>${workflow.runName}</td></tr>
                <tr><td>Status</td><td class="${workflow.success ? 'success' : 'failed'}">${status}</td></tr>
                <tr><td>Completed</td><td>${workflow.complete}</td></tr>
                <tr><td>Duration</td><td>${workflow.duration}</td></tr>
                <tr><td>Exit Status</td><td>${workflow.exitStatus}</td></tr>
                <tr><td>Error Report</td><td>${workflow.errorReport ?: 'None'}</td></tr>
                <tr><td>Command</td><td><code>${workflow.commandLine}</code></td></tr>
                <tr><td>Project Dir</td><td>${workflow.projectDir}</td></tr>
                <tr><td>Work Dir</td><td>${workflow.workDir}</td></tr>
                <tr><td>Output Dir</td><td>${params.outdir}</td></tr>
            </table>
            
            <h3>Parameters Used</h3>
            <pre>${summary_log}</pre>
            
            <p>
                <small>
                    This email was automatically generated by Nextflow<br>
                    ${workflow.manifest.name} v${workflow.manifest.version}<br>
                    ${workflow.manifest.homePage}
                </small>
            </p>
        </body>
        </html>
        """
        
        return html
    }
}