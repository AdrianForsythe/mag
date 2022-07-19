process taxa_list {
    label 'process_low'
    publishDir path: "$params.outdir/taxa_list", mode: publish_dir_mode

    input:
    path(taxa)

    output:
    file("*.decontam-children.fastq.gz"), emit: taxa_list
    path "versions.yml"                                   , emit: versions

    when:
    params.decontamiante

    script:
    out = "${taxa.basename}.tsv"
    """
    cat $taxa | tr "," "\t" | awk '{print \$1}' | grep -v "Tax_ID" | tr "\n" " " > $out
    """
}
