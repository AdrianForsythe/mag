process decontaminate {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::krakentools=1.2" : null)
    publishDir path: "$params.outdir/decontaminate", mode: publish_dir_mode

    input:
    tuple val(meta.id), val(meta), path(fastq), path(k_report)
    path(taxa)

    output:
    tuple val(meta), file("*.decontam-children.fastq.gz"), emit: decontam_output
    path "versions.yml"                                   , emit: versions

    when:
    params.decontamiante

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_kraken_reads.py -k ${k_report} -s ${fastq} -o ${prefix} -t ${taxa} --fastq-output --exclude
    pigz -f -p ${tasks.cpu} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        KrakenTools: KrakenTools version 1.2 Copyright 2021 Jennifer Lu
    END_VERSIONS

    """
}
