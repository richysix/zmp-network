def convert_species_lc_nospace (species) {
    species_nospace = species.replaceAll(" ", "_")
    species_lc = species_nospace.toLowerCase()
    return species_lc
}

process GET_ANNOTATION {
    // tag "$meta.id"
    label 'local'
    publishDir "reference", pattern: "$outFile"

    input:
    val(species)
    val(ensemblVersion)

    output:
    path("$outFile")

    script:
    species_lc_nospace = convert_species_lc_nospace(species)
    outFile = species_lc_nospace + "-e" + ensemblVersion + "-annotation.txt"
    // convert species to lower case and remove spaces
    """
    # Download get annotation scripts
    wget $params.getAnnoScriptURL
    wget $params.getAnnoBashURL
    # make it executable
    chmod a+x get_ensembl_gene_annotation.pl get_ensembl_gene_annotation.sh
    # Run 
    ./get_ensembl_gene_annotation.sh -s "$species" -e $ensemblVersion -o $outFile
    """
}

process GET_GO_ANNOTATION {
    // tag "$meta.id"
    label 'local'
    publishDir "reference", pattern: "*_e*_go.txt"

    input:
    val(species)
    val(ensemblVersion)

    output:
    path("*_e*_go.txt")

    script:
    species_lc_nospace = convert_species_lc_nospace(species)
    """
    # Download get annotation scripts
    wget $params.getGOAnnoBashURL
    # make it executable
    chmod a+x get-go-annotation.sh
    ./get-go-annotation.sh -s $species_lc_nospace -e $ensemblVersion
    """
}

// Workflow for testing the module
// nextflow.enable.moduleBinaries = true
params.getGOAnnoBashURL="https://github.com/richysix/uge-job-scripts/raw/e5f5faad28b23ff55419726beef3675f9e5fdba3/get-go-annotation.sh"
params.EnsemblVersionGO="105"
workflow GET_GO_ANNOTATION_WF {
    GET_GO_ANNOTATION(params.Species, params.EnsemblVersionGO)
        .view()
}

// Workflow for testing the module
params.getAnnoScriptURL="https://github.com/iansealy/bio-misc/raw/4e27d60323907d55d37a1ec8f468ca771f542f78/get_ensembl_gene_annotation.pl"
params.getAnnoBashURL="https://github.com/richysix/uge-job-scripts/raw/e5f5faad28b23ff55419726beef3675f9e5fdba3/get_ensembl_gene_annotation.sh"
params.EnsemblVersion="100"
params.Species='Mus musculus'
workflow GET_ANNOTATION_WF {
    GET_ANNOTATION(params.Species, params.EnsemblVersion)
        .view()
}

// For testing
workflow {
    GET_ANNOTATION_WF()

    GET_GO_ANNOTATION_WF()
}