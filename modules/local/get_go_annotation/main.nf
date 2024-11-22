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
    """
    ${params.ScriptDir}/get-go-annotation.sh -s $species -e $ensemblVersion
    """
}

// Workflow for testing the module
// nextflow.enable.moduleBinaries = true
params.ScriptDir="$HOME/checkouts/zmp-network/bin"
params.EnsemblVersionGO=100
params.Species='mus_musculus'
workflow {
    GET_GO_ANNOTATION(params.Species, params.EnsemblVersionGO)
        .view()
}