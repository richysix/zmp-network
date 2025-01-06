def convert_species_lc_nospace (species) {
    def species_nospace = species.replaceAll(" ", "_")
    def species_lc = species_nospace.toLowerCase()
    return species_lc
}

process GET_ANNOTATION {
    // tag "$meta.id"
    label 'process_single'

    input:
    val(species)
    val(ensemblVersion)

    output:
    path("$out_file")

    script:
    // convert species to lower case and remove spaces
    species_lc_nospace = convert_species_lc_nospace(species)
    out_file = species_lc_nospace + "-e" + ensemblVersion + "-annotation.txt"
    """
    module load Ensembl/${ensemblVersion}
    get_ensembl_gene_annotation.pl --species "${species}" > ${out_file}
    """
}

// Workflow for testing the GET_ANNOTATION process
// params.ensembl_version="100"
// params.species='danio_rerio'
workflow GET_ANNOTATION_WF {
    GET_ANNOTATION(
        params.species,
        params.ensembl_version
    ).view()
}
