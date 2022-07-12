process CENTRIFUGE_DB_PREPARATION {

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path db

    output:
    tuple val("${db.toString().replace(".tar.gz", "")}"), path("*.cf"), emit: db
    path "versions.yml"                                               , emit: versions

    when:
    !params.skip_taxonomic_classification

    script:
    """
    tar -xf "${db}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n 1p | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
