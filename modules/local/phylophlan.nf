process PHYLOPHLAN {
    tag "$species"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::phylophlan, bioconda::diamond, bioconda::blast, bioconda::FastTree" : null)

    input:
    tuple val(species), val(genus), val(outgroup)
    path(bins).splitCsv(header: ['user_genome','classification','fastani_reference','fastani_reference_radius','fastani_taxonomy','fastani_ani','fastani_af','closest_placement_reference','closest_placement_taxonomy','closest_placement_ani','closest_placement_af','pplacer_taxonomy','classification_method','note','other_related_references','msa_percent','red_value','warnings'], skip: 1 ).filter("$species")[0]
    tuple val(bin_name).filter("$bins")
    path(fnas)

    output:
    tuple val(genus), val(species), path("${prefix}/*.fa"), emit: fa
    path "phylophlan/$species/tree/*.tre", emit: tree
    path "phylophlan/$species/input/*.[fa,fna]", emit: fastas
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    nref = params.phylophlan_nref
    diversity = params.phylophlan_diversity
    outdir_species = "phylophlan/$species"
    species_database = "$outdir_species/database"
    indir = "$outdir_species/input/"
    outdir_tree = "$outdir_species/tree/"
    outdir_outgroup = "$outdir_species/outgroup/"

    // def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    // def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    """
    mkdir "$outdir_species"
    mkdir "$indir"
    mkdir "$outdir_tree"
    mkdir "$outdir_outgroup"

    # make links to input bins
    for i in $fnas
    do
        prefix = basename $i .fna
        ln -s -f $i ${indir}/${prefix}.fna
    done

    # Automatically retrieve the set of core UniRef90 proteins for each species:
    phylophlan_setup_database \
        -g $species \
        -o $outdir_species \
        --overwrite \
        --verbose 2>&1 | tee $outdir_species/phylophlan_setup_database.log

    # If we want to add public available reference genomes, we can use PhyloPhlAn 3.0 to download them:
    phylophlan_get_reference \
        -g $species \
        -o $species_database \
        --database_update \
        -m ${outdir_species}/assembly_summary_genbank.txt \
        -n $nref \
        --verbose 2>&1 | tee $outdir_species/phylophlan_get_reference.log

    # do the same for the outgroup
    phylophlan_get_reference \
        -g $outgroup \
        -o ${outdir_outgroup}/ \
        --database_update \
        -m ${outdir_species}/assembly_summary_genbank.txt \
        -n 1 \
        --verbose 2>&1 | tee $outdir_species/phylophlan_get_outgroup.log

    # copy outgroup to input folder
    outgroup_NAME=$(find $outdir_outgroup/* |  awk '{print \$1}'|  sed 's!.*/!!'| sed 's/.fna.gz//g')
    pigz -f -p $task.cpus $outdir_outgroup/*.fna.gz
    cp $outdir_outgroup/*.fna $indir/
    pigz -f -p $task.cpus $outdir_outgroup/*.fna

    # configuration file generated using the following command:
    phylophlan_write_config_file \
        -o "${outdir_species}"/references_config.cfg \
        -d a \
        --db_aa diamond \
        --map_aa diamond \
        --map_dna diamond \
        --msa mafft \
        --overwrite \
        --trim trimal \
        --tree1 fasttree

    # and we can now build a phylogenetic tree of all our genomes:
    # make sure to specify where the database file ends up!

    phylophlan \
        -i ${indir}/ \
        --output_folder ${outdir_tree} \
        -d $species_database \
        -t a \
        -f ${outdir_species}/references_config.cfg \
        --nproc $task.cpus \
        --diversity $diversity \
        --fast \
        --verbose 2>&1 | tee ${outdir_species}/phylophlan_${species}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phylophlan: \$(echo \$(phylophlan --version 2>&1) | sed 's/^.*phylophlan //')
    END_VERSIONS
    """
}
