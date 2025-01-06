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

