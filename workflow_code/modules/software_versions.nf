process SOFTWARE_VERSIONS {
    input:
        path(versions_file)
    
    output:
        path "software_versions${params.assay_suffix}.md", emit: software_versions
        path "software_versions${params.assay_suffix}.yaml", emit: software_versions_yaml

    script:
    """
    software_versions.py ${versions_file} software_versions${params.assay_suffix}.md --workflow NF_RCP --workflow_version ${workflow.manifest.version} --assay rnaseq
    """
}