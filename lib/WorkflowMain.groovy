//
// This file holds several functions specific to the main.nf workflow in the h3abionet/chipimputation pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The H3ABioNet Consortium\n" +
            "  https://doi.org/10.1038/s41467-018-05188-3\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //
    public static String toolCitationText(params) {

        def citation_text = [
            "Tools used in the workflow included:",
            "Eagle (Loh et al., 2016),",
            "Minimac4 (Das et al., 2016),",
            "BCFtools (Danecek et al., 2021),",
            "PLINK (Purcell et al., 2007).",
            "Pipeline was executed with Nextflow version ${workflow.nextflow.version} (Di Tommaso et al., 2017)",
            "with wrappers available from",
            "https://github.com/h3abionet/chipimputation."
        ].join(' ').trim()

        return citation_text
    }

    //
    // Generate workflow version string
    //
    public static String getWorkflowVersion() {
        String version_string = ""
        if (workflow.manifest.version) {
            def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
            version_string += "${prefix_v}${workflow.manifest.version}"
        }

        if (workflow.commitId) {
            def git_shortsha = workflow.commitId.substring(0, 7)
            version_string += "-g${git_shortsha}"
        }

        return version_string
    }

    //
    // Get software versions
    //
    public static String processVersionsFromYAML(yaml_file) {
        Yaml yaml = new Yaml()
        versions = yaml.load(yaml_file).collectEntries { k, v -> [ k.tokenize(':')[-1], v ] }
        return yaml.dumpAsMap(versions).trim()
    }
}