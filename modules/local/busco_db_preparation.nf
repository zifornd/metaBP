process BUSCO_DB_PREPARATION {
    tag "${database.baseName}"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"
    input:
    path database

    output:
    path "buscodb/*", emit: db
    path database

    script:
    """
    mkdir buscodb
    tar -xf ${database} -C buscodb
    """
}
