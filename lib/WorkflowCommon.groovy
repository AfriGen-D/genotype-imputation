class WorkflowCommon {

    //
    // Get genome attribute
    //
    public static getGenomeAttribute(params, attribute) {
        def val = null
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[params.genome].containsKey(attribute)) {
                val = params.genomes[params.genome][attribute]
            }
        }
        return val
    }

    //
    // Get workflow name
    //
    public static getWorkflowName(workflow, manifest) {
        def workflow_name = manifest.name ?: workflow.manifest.name
        return workflow_name
    }

    //
    // Get workflow version
    //
    public static getWorkflowVersion(workflow, manifest) {
        def workflow_version = manifest.version ?: workflow.manifest.version
        return workflow_version
    }

    //
    // Get workflow description
    //
    public static getWorkflowDescription(workflow, manifest) {
        def workflow_description = manifest.description ?: workflow.manifest.description
        return workflow_description
    }

    //
    // Check if a parameter is defined
    //
    public static hasParam(params, param_name) {
        return params.containsKey(param_name) && params[param_name]
    }

    //
    // Create parameter summary for logging
    //
    public static paramsSummaryMap(workflow, params, log) {
        def summary = [:]
        def header = "----------------------------------------------------"
        
        summary['Pipeline Name'] = workflow.manifest.name
        summary['Pipeline Version'] = workflow.manifest.version
        summary['Run Name'] = workflow.runName
        
        // Core parameters
        if (params.input) summary['Input'] = params.input
        if (params.outdir) summary['Output dir'] = params.outdir
        if (params.genome) summary['Genome'] = params.genome
        
        // Resources
        summary['Max Resources'] = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"
        
        // Container
        if (workflow.containerEngine) summary['Container'] = "${workflow.containerEngine} - ${workflow.container}"
        
        // Profile
        summary['Profile'] = workflow.profile
        
        return summary
    }

    //
    // Check and return file if exists
    //
    public static returnFile(it) {
        if (!file(it).exists()) exit 1, "File not found: ${it}"
        return file(it)
    }

    //
    // Function to check max resources
    //
    public static checkMaxResources(params, log) {
        // Memory
        if (params.max_memory) {
            def max_memory = params.max_memory as nextflow.util.MemoryUnit
            if (max_memory > 128.GB) {
                log.warn "Max memory set to > 128GB. Ensure your system can handle this."
            }
        }
        
        // CPUs
        if (params.max_cpus > 64) {
            log.warn "Max CPUs set to > 64. Ensure your system can handle this."
        }
        
        // Time
        if (params.max_time) {
            def max_time = params.max_time as nextflow.util.Duration
            if (max_time > 240.h) {
                log.warn "Max time set to > 240 hours. Jobs may timeout on some systems."
            }
        }
    }

    //
    // Check if running offline
    //
    public static isOffline(workflow, params) {
        try {
            return workflow.offline || params.offline
        } catch(Exception e) {
            return false
        }
    }

    //
    // Check if container is available
    //
    public static checkContainer(workflow, params, log) {
        if (workflow.profile.contains('docker') || workflow.profile.contains('singularity')) {
            if (!workflow.containerEngine) {
                log.error "Container engine not detected. Please install Docker or Singularity."
                System.exit(1)
            }
        }
    }

    //
    // Get file extension
    //
    public static getFileExtension(filename) {
        def extension = ''
        def i = filename.lastIndexOf('.')
        if (i > 0) {
            extension = filename.substring(i+1)
        }
        return extension
    }

    //
    // Check mandatory parameters
    //
    public static checkMandatoryParams(params, mandatoryParams, log) {
        def missingParams = []
        mandatoryParams.each { param ->
            if (!params.containsKey(param) || !params[param]) {
                missingParams.add(param)
            }
        }
        
        if (missingParams.size() > 0) {
            log.error "The following mandatory parameters are missing: ${missingParams.join(', ')}"
            System.exit(1)
        }
    }
}