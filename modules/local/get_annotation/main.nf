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

process GET_GO_ANNOTATION {
    // tag "$meta.id"
    label 'process_single'
    publishDir "reference", pattern: "*_e*_go.txt"

    input:
    val(species)
    val(ensemblVersionGO)

    output:
    path("*_e*_go.txt")

    script:
    species_lc_nospace = convert_species_lc_nospace(species)
    """
    get-go-annotation.sh -s $species_lc_nospace -e $ensemblVersionGO
    """
}

// Workflow for testing the GET_GO_ANNOTATION process
// params.ensembl_versionGO="105"
workflow GET_GO_ANNOTATION_WF {
    GET_GO_ANNOTATION(
        params.species,
        params.ensembl_versionGO
    ).view()
}

// For testing
workflow {
    GET_ANNOTATION_WF()

    GET_GO_ANNOTATION_WF()
}