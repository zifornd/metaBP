// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process GTDBTK_DB_PREPARATION {
    tag "${database}"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ubuntu:20.04"
    } else {
        container "ubuntu:20.04"
    }

    input:
    path(database)

    output:
    tuple val("${database.toString().replace(".tar.gz", "")}"), path("database/*")

    script:
    """
    mkdir database
    tar -xzf ${database} -C database --strip 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n 1p | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
