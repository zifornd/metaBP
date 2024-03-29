process KRAKEN2_DB_PREPARATION {
    
    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        "biocontainers/biocontainers:v1.2.0_cv1" }"

    input:
    path db

    output:
    tuple val("kraken2_db"), path("database/*.k2d"), emit: db

    script:
    """
    mkdir db_tmp
    tar -xf "${db}" -C db_tmp
    mkdir database
    mv `find db_tmp/ -name "*.k2d"` database/
    """
}
