include { convert_species_lc_nospace } from '../../../modules/local/get_annotation/main'

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

