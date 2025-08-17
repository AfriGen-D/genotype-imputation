//
// This file holds several utility functions used within the pipeline
//

import org.yaml.snakeyaml.Yaml
import groovy.json.JsonSlurper

class Utils {

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by priority: first items have higher priority
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that channels are in the right order
        def channels_in_wrong_order = false
        required_channels_in_order.each { channel ->
            if (channels.indexOf(channel) > required_channels_in_order.indexOf(channel)) {
                channels_in_wrong_order = true
            }
        }

        if (channels_missing | channels_in_wrong_order) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

    //
    // Function to extract metadata from a sample sheet
    //
    public static LinkedHashMap parseSampleSheet(sample_sheet_path) {
        def sample_sheet = sample_sheet_path.toString()
        def input_rows = []
        def header = []
        
        // Read the CSV file
        new File(sample_sheet).withReader { reader ->
            def lines = reader.readLines()
            header = lines[0].tokenize(',')
            
            lines[1..-1].each { line ->
                def row = line.tokenize(',')
                def sample_map = [:]
                header.eachWithIndex { col, idx ->
                    sample_map[col] = row[idx]
                }
                input_rows.add(sample_map)
            }
        }
        
        return input_rows
    }

    //
    // Function to check if a file exists
    //
    public static Boolean fileExists(file_path) {
        def file = new File(file_path.toString())
        return file.exists() && file.isFile()
    }

    //
    // Function to check if all required params are provided
    //
    public static void checkParams(params, log) {
        def required = ['input', 'outdir']
        def missing = []
        
        required.each { param ->
            if (!params[param]) {
                missing.add(param)
            }
        }
        
        if (missing.size() > 0) {
            log.error "The following required parameters are missing: ${missing.join(', ')}"
            System.exit(1)
        }
    }

    //
    // Function to initialise default params
    //
    public static Map initialiseParams(params, log) {
        def new_params = params.clone()
        
        // Set defaults if not provided
        if (!new_params.chromosomes) new_params.chromosomes = 'ALL'
        if (!new_params.chunk_size) new_params.chunk_size = 5000000
        if (!new_params.minRatio) new_params.minRatio = 0.01
        if (!new_params.eagle_pbwt_iters) new_params.eagle_pbwt_iters = 2
        if (!new_params.min_ac) new_params.min_ac = 2
        if (!new_params.site_missingness) new_params.site_missingness = 0.05
        
        return new_params
    }

    //
    // Function to save parameters to JSON for reporting
    //
    public static void saveParams(params, output_dir) {
        def json_str = new groovy.json.JsonBuilder(params).toPrettyString()
        def json_file = new File("${output_dir}/pipeline_info/params.json")
        json_file.parentFile.mkdirs()
        json_file.text = json_str
    }
}