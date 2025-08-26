/*
========================================================================================
    COMMON WORKFLOW FUNCTIONS
========================================================================================
*/

def getGenomeAttribute(params, attribute) {
    def val = null
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            val = params.genomes[params.genome][attribute]
        }
    }
    return val
}

def getWorkflowName(workflow, manifest) {
    def workflow_name = manifest.name ?: workflow.manifest.name
    return workflow_name
}

def getWorkflowVersion(workflow, manifest) {
    def workflow_version = manifest.version ?: workflow.manifest.version
    return workflow_version
}

def checkMaxResources(params, log) {
    // Memory
    if (params.max_memory) {
        // Just check if it's a large value string
        if (params.max_memory.toString().contains('256') || params.max_memory.toString().contains('512')) {
            log.warn "Max memory set to very high value. Ensure your system can handle this."
        }
    }
    
    // CPUs
    if (params.max_cpus && params.max_cpus > 64) {
        log.warn "Max CPUs set to > 64. Ensure your system can handle this."
    }
    
    // Time
    if (params.max_time) {
        // Just check if it contains large hour values
        if (params.max_time.toString().contains('240') || params.max_time.toString().contains('480')) {
            log.warn "Max time set to very long duration. Jobs may timeout on some systems."
        }
    }
}

def isOffline(workflow, params) {
    try {
        return workflow.offline || params.offline
    } catch(Exception e) {
        return false
    }
}

def checkContainer(workflow, params, log) {
    if (workflow.profile.contains('docker') || workflow.profile.contains('singularity')) {
        if (!workflow.containerEngine) {
            log.warn "Container engine not detected but profile requires it."
        }
    }
}

def getFileExtension(filename) {
    def extension = ''
    def i = filename.lastIndexOf('.')
    if (i > 0) {
        extension = filename.substring(i+1)
    }
    return extension
}

def checkMandatoryParams(params, mandatoryParams, log) {
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

def createWorkflowSummary(workflow, params) {
    def summary = [:]
    summary['version'] = workflow.manifest.version
    summary['runName'] = workflow.runName
    summary['success'] = workflow.success
    summary['dateComplete'] = workflow.complete
    summary['duration'] = workflow.duration
    summary['exitStatus'] = workflow.exitStatus
    summary['errorMessage'] = (workflow.errorMessage ?: 'None')
    summary['errorReport'] = (workflow.errorReport ?: 'None')
    summary['commandLine'] = workflow.commandLine
    summary['projectDir'] = workflow.projectDir
    summary['workDir'] = workflow.workDir
    summary['launchDir'] = workflow.launchDir
    summary['configFiles'] = workflow.configFiles.join(', ')
    summary['profile'] = workflow.profile
    summary['container'] = workflow.container
    summary['containerEngine'] = workflow.containerEngine
    summary['resume'] = workflow.resume
    
    return summary
}