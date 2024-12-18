//
// Subworkflow to run GET_ANNOTATION and GET_GO_ANNOTATION
//

include { GET_ANNOTATION     } from '../../modules/local/get_annotation/main'
include { GET_GO_ANNOTATION  } from '../../modules/local/get_annotation/main'
include { convert_species_lc_nospace } from '../../modules/local/get_annotation/main'

workflow GET_ANNO_GET_GO_ANNO {
    take:
    ch_species    // channel: species name from params.species
    ch_anno       // channel: [ anno_file, version, anno_script_url, anno_bash_url ]
    ch_go_anno    // channel: [ go_anno_file, go_version, go_bash_url ]
    
    main:
    // Check if annotation file already exists
    if (ch_anno.anno_file.exists()) {
        ch_anno_file = Channel.of(ch_anno.anno_file)
        if ( params.debug > 1) {
            println("Annotation file exists")
            ch_anno_file.view { x -> "Annotation file: $x" }
        }
    } else {
        ch_anno_file = GET_ANNOTATION(
            ch_species,
            ch_anno.version,
            ch_anno.anno_script_url,
            ch_anno.anno_bash_url
        )
        if ( params.debug > 1) {
            println("Annotation file does NOT exist. Downloading...")
            ch_anno_file.view { x -> "Annotation file: $x" }
        }
    }

    if (ch_go_anno.go_anno_file.exists()) {
        ch_go_anno_file = Channel.of(ch_go_anno.go_anno_file)
        if ( params.debug > 1) {
            println("GO annotation file exists")
            ch_go_anno_file.view { x -> "GO annotation file: $x" }
        }
    } else {
        ch_go_anno_file = GET_GO_ANNOTATION(
            ch_species,
            ch_go_anno.go_version,
            ch_go_anno.go_bash_url
        )
        if ( params.debug > 1) {
            println("GO annotation file does NOT exist. Downloading...")
            ch_go_anno_file.view { x -> "GO annotation file: $x" }
        }
    }

    emit:
    // convert to value channels using .first()
    anno_file = ch_anno_file        // value channel: [ path(annotation_file) ]
    go_anno_file = ch_go_anno_file  // value channel: [ path(GO_annotation_file) ]
}

params.get_anno_script_url="https://github.com/iansealy/bio-misc/raw/4e27d60323907d55d37a1ec8f468ca771f542f78/get_ensembl_gene_annotation.pl"
params.get_anno_bash_url="https://github.com/richysix/uge-job-scripts/raw/e5f5faad28b23ff55419726beef3675f9e5fdba3/get_ensembl_gene_annotation.sh"
params.ensembl_version="100"
params.species='danio_rerio'
params.get_go_anno_bash_url="https://github.com/richysix/uge-job-scripts/raw/e5f5faad28b23ff55419726beef3675f9e5fdba3/get-go-annotation.sh"
params.ensembl_versionGO="105"
workflow {
    GET_ANNO_GET_GO_ANNO(
        params.species,
        [ anno_file: file("reference/danio_rerio-e100-annotation.txt"), 
          version: params.ensembl_version,
          anno_script_url: params.get_anno_script_url,
          anno_bash_url: params.get_anno_bash_url ],
        [ go_anno_file: file("reference/danio_rerio_e105_go.txt"),
          go_version: params.ensembl_versionGO,
          go_bash_url: params.get_go_anno_bash_url ]
    )
    
    GET_ANNO_GET_GO_ANNO.out.anno_file.view()
    GET_ANNO_GET_GO_ANNO.out.go_anno_file.view()
}