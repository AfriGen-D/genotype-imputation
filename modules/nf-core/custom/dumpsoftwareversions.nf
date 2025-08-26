process CUSTOM_DUMPSOFTWAREVERSIONS {
    tag "versions"
    label 'process_single'
    
    input:
    path versions
    
    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions.tsv"    , emit: tsv
    path "versions.yml"            , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    #!/usr/bin/env python3
    import yaml
    import platform
    from textwrap import dedent
    
    def _make_versions_html(versions):
        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions">
                    <thead>
                        <tr>
                            <th> Process Name </th>
                            <th> Software </th>
                            <th> Version  </th>
                        </tr>
                    </thead>
                '''
            )
        ]
        for process, process_versions in versions.items():
            html.append("<tbody>")
            for i, (tool, version) in enumerate(process_versions.items()):
                html.append(
                    dedent(
                        f'''\\
                        <tr>
                            <td>{process if i == 0 else ""}</td>
                            <td>{tool}</td>
                            <td>{version}</td>
                        </tr>
                        '''
                    )
                )
            html.append("</tbody>")
        html.append("</table>")
        return "\\n".join(html)
    
    versions_this_module = {}
    versions_this_module["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    
    with open("$versions") as f:
        versions_by_process = yaml.load(f, Loader=yaml.BaseLoader) | versions_this_module
    
    # aggregate versions by the module name (derived from fully-qualified process name)
    versions_by_module = {}
    for process, process_versions in versions_by_process.items():
        module = process.split(":")[-1]
        try:
            assert versions_by_module[module] == process_versions, (
                "We assume that software versions are the same between all modules. "
                "If you see this error-message it means you discovered an edge-case "
                "and should open an issue on GitHub. "
            )
        except KeyError:
            versions_by_module[module] = process_versions
    
    versions_by_module["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "$workflow.manifest.name": "$workflow.manifest.version",
    }
    
    versions_mqc = {
        "id": "software_versions",
        "section_name": "${workflow.manifest.name} Software Versions",
        "section_href": "https://github.com/${workflow.manifest.name}",
        "plot_type": "html",
        "description": "are collected at run time from the software output.",
        "data": _make_versions_html(versions_by_module),
    }
    
    with open("software_versions.yml", "w") as f:
        yaml.dump(versions_by_module, f, default_flow_style=False)
    with open("versions.yml", "w") as f:
        yaml.dump(versions_this_module, f, default_flow_style=False)
    
    # Create TSV file
    with open("software_versions.tsv", "w") as f:
        f.write("Process\\tSoftware\\tVersion\\n")
        for process, tools in versions_by_module.items():
            for tool, version in tools.items():
                f.write(f"{process}\\t{tool}\\t{version}\\n")
    """
    
    stub:
    """
    touch software_versions.yml
    touch software_versions.tsv
    touch versions.yml
    """
}