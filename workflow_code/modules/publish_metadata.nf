process PUBLISH_METADATA {
    input:
        path(runsheet)
        path(isa_archive)

    output:
        path("metadata/*")

    script:
        def copy_isa = isa_archive ? "cp -L '${isa_archive}' metadata/" : ''
        """
        mkdir -p metadata
        cp -L '${runsheet}' metadata/
        ${copy_isa}
        """
}
