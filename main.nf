#!/usr/bin/env nextflow

log.info """\
  NETWORK CONSTRUCTION PIPELINE
  -----------------------------

  Testing: ${params.testing}
  Clustering: ${params.clustering}
  Debug: ${params.debug}
  Expt to sample file: ${params.expts}
  All counts file: ${params.all_counts}
  Threshold params: ${params.threshold}
  KNN params: ${params.knn}
  Inflation Params: ${params.inflation_params}
"""

include { GET_ANNO_GET_GO_ANNO } from './subworkflows/local/get_gene_and_go_annotation'

def get_threshold(m) {
    def t
    if (m[0][1] == "t") {
        t = (m[0][2].toInteger()) / 100
    } else {
        t = 0
    }
    return t
}

// process to subset the file containing all the samples
// to each separate experiment
// see README
process SUBSET {
    label 'big_mem_retry'
    publishDir "results", pattern: "expts.txt"

    input: 
    val expt_file
    val all_counts_file

    output:
    tuple(
        path("expts.txt"),
        path("*/samples.tsv"),
        path("*/counts-by-gene.tsv")
    )

    script:
    debug = params.debug ? "-d" : ""
    """
    subset-by-expt.sh ${debug} ${expt_file} ${all_counts_file}
    """
}

process CREATE_BASE_NETWORK {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/*all-tpm*"

    input: 
    tuple val(expt_dir), path(sample_file), path(count_file)

    output:
    tuple val(expt_dir),
        path("${expt_dir}/all-tpm.tsv"), path("${expt_dir}/all-tpm.tab"),
        path("${expt_dir}/all-tpm-orig.mcx"), path("${expt_dir}/all-tpm-t20.mcx"),
        path("${expt_dir}/all-tpm-orig.mat.csv.gz"),
        path("${expt_dir}/${expt_dir}-all-tpm-orig.cor-hist.txt")

    script:
    """
    mkdir -p ${expt_dir}
    awk -F"\\t" '{if(NR > 1){ print \$2 "\\t" \$3 }}' \
     ${sample_file} > ${expt_dir}/samples.txt

    module load R/${params.r_version}
    counts-to-fpkm-tpm.R \
    --transcripts_file $params.transcript_file \
    --output_base ${expt_dir}/all --output_format tsv \
    --tpm ${expt_dir}/samples.txt $count_file

    module load MCL/$params.mcl_version

    # make network with all edges in
    mcxarray -data ${expt_dir}/all-tpm.tsv -co 0 \
    $params.skip_rows $params.skip_cols \
    $params.cor_measure $params.labels \
    --write-binary -o ${expt_dir}/all-tpm-orig.mcx \
    -write-tab ${expt_dir}/all-tpm.tab

    mcxdump -imx ${expt_dir}/all-tpm-orig.mcx \
    -tab ${expt_dir}/all-tpm.tab --dump-table \
    -digits 3 -sep-field "," -sep-lead "," \
    -o ${expt_dir}/all-tpm-orig.mat.csv
    gzip ${expt_dir}/all-tpm-orig.mat.csv

    module load Python/$params.python_version
    cor-hist.py ${expt_dir}/all-tpm-orig.mat.csv.gz \
    ${expt_dir}/${expt_dir}-all-tpm-orig.cor-hist.txt

    # Also make one filtered with abs(), gt > 0.2
    mcx alter -imx ${expt_dir}/all-tpm-orig.mcx \
    -tf "abs(), gt(0.2)" \
    --write-binary -o ${expt_dir}/all-tpm-t20.mcx
    """
}

// Create a basic network with a correlation threshold of 0.2
// Test varying threshold and knn parameters
process TEST_PARAMETERS {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(mci_file)

    output:
    tuple val(dir), path("$dir/$dir-all-tpm.cor-stats.tsv"),
        path("$dir/$dir-all-tpm.knn-stats.tsv")

    script:
    """
    module load MCL/$params.mcl_version
    
    mkdir -p $dir

    # vary correlation
    mcx query -imx ${mci_file} --vary-correlation \
    --output-table > $dir/$dir-all-tpm.cor-stats.tsv

    # test varying k-nearest neighbours
    mcx query -imx ${mci_file} -vary-knn $params.knn_test_params \
    --output-table > $dir/$dir-all-tpm.knn-stats.tsv
    """
}

// Threshold
process FILTER_COR {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(mci_file)
    each threshold

    output:
    tuple val(dir), path("$dir/*-[t][0-9]*.mcx"),
        path("$dir/$dir-all-tpm-t*.stats.tsv")

    script:
    Integer suffix = threshold * 100
    """
    module load MCL/$params.mcl_version

    mkdir -p $dir
    mcx alter -imx ${mci_file} \
    -tf "gt(${threshold}), add(-${threshold})" \
    --write-binary -o $dir/all-tpm-t${suffix}.mcx
    mcx alter -imx $dir/all-tpm-t${suffix}.mcx \
    -tf "add(${threshold})" | mcx query -imx - > \
    $dir/$dir-all-tpm-t${suffix}.stats.tsv
    """
}

process FILTER_KNN {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"

    input:
    tuple val(dir), path(mci_file)
    each knn_threshold

    output:
    tuple val(dir), path("$dir/*-[k][0-9]*.mcx"),
        path("$dir/$dir-all-tpm-t*-k*.stats.tsv")

    script:
    threshold=0.2
    thresholdSuffix=20
    """
    module load MCL/$params.mcl_version

    mkdir -p $dir
    knnBase="all-tpm-t${thresholdSuffix}-k${knn_threshold}"
    mcx alter -imx ${mci_file} \
    -tf "abs(), gt(${threshold}), add(-${threshold}), #knn($knn_threshold)" \
    --write-binary -o $dir/\${knnBase}.mcx
    # stats
    mcx alter -imx $dir/\${knnBase}.mcx -tf "add(${threshold})" | \
    mcx query -imx - > $dir/$dir-\${knnBase}.stats.tsv
    """
}

process FILTER_STATS {
    label 'retry'
    publishDir "results", pattern: "plots/*.pdf"

    input:
    path("expts.txt")
    path("*")
    path("*")
    path("*")

    output:
    path("plots/*.pdf")

    script:
    """
    module load R/${params.r_version}
    edge-filtering-analysis.R --samples_file $params.expts
    """
}

// Cluster network
process CLUSTER {
    label 'med_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(mci_file), path(stats_file)
    each inflation

    output:
    tuple val(dir), path(mci_file), path("*/*.mcx.I[0-9][0-9]"),
        path("*/*.mcx.I[0-9][0-9].stats.tsv"), path("*/*.mcx.I[0-9][0-9].cl*")

    script:
    Integer inflationSuffix = inflation * 10
    """
    mkdir -p ${dir}

    module load MCL/${params.mcl_version}
    mcl $mci_file -I $inflation -o $dir/${mci_file}.I${inflationSuffix}
    
    mci_base=\$( basename $mci_file .mci )
    clm info $mci_file $dir/${mci_file}.I${inflationSuffix} >> $dir/\${mci_base}.info.txt

    clm info --node-all-measures --node-self-measures $mci_file \
    $dir/${mci_file}.I${inflationSuffix} > $dir/${mci_file}.I${inflationSuffix}.stats.tsv

    module load Python/$params.python_version
    summarise_clustering.py --expt_name $dir \
    $dir/${mci_file}.I${inflationSuffix}
    """
}

process GBA {
    label 'huge_mem_retry'
    publishDir "results", pattern: "*/all-tpm*.{graphml,tsv,csv,pdf}"

    input:
    tuple val(dir), path(tab_file), path(mci_file), path(cluster_file)
    path(annotation_file)
    path(go_annotation_file)
    path(zfa_annotation_file)

    output:
    tuple val(dir), path(cluster_file), path("*/all-tpm*.graphml"),
        path("*/all-tpm*.nodes.tsv"), path("*/all-tpm*.edges.tsv"),
        path("*/$dir-all-tpm*.auc.tsv"), path("*/all-tpm*.gene-scores.tsv"),
        path("*/all-tpm*.GBA-plots.pdf")

    script:
    matches = (mci_file =~ /all-tpm-(t?)(\d*)-?(k?)(\d*).mcx$/)
    t_num = get_threshold(matches)
    if (params.debug){
        println("GBA: Threshold value = " + t_num)
    }
    """
    mkdir -p ${dir}
    cluster_base=\$( basename $cluster_file )

    # convert from binary
    module load MCL/${params.mcl_version}
    mci_base=\$( basename $mci_file .mcx)
    mcx convert $mci_file \${mci_base}.mci

    # run convert_mcl script
    module load Python/$params.python_version
    convert_mcl.py \
    --min_cluster_size 4 --graph_id \${cluster_base} \
    --graphml_file ${dir}/\${cluster_base}.graphml \
    --nodes_file ${dir}/\${cluster_base}.nodes.tsv \
    --edges_file ${dir}/\${cluster_base}.edges.tsv \
    --edge_offset ${t_num} \
    \${mci_base}.mci $cluster_file $tab_file $annotation_file

    module load R/${params.r_version}
    run-GBA-network.R \
    --auc_file ${dir}/${dir}-\${cluster_base}.go.auc.tsv \
    --scores_file ${dir}/\${cluster_base}.go.gene-scores.tsv \
    --plots_file ${dir}/\${cluster_base}.go.GBA-plots.pdf \
    --min.term.size $params.min_term_size --max.term.size $params.max_term_size \
    ${dir}/\${cluster_base}.nodes.tsv \
    ${dir}/\${cluster_base}.edges.tsv \
    $go_annotation_file

    run-GBA-network.R \
    --auc_file ${dir}/${dir}-\${cluster_base}.zfa.auc.tsv \
    --scores_file ${dir}/\${cluster_base}.zfa.gene-scores.tsv \
    --plots_file ${dir}/\${cluster_base}.zfa.GBA-plots.pdf \
    --min.term.size $params.min_term_size --max.term.size $params.max_term_size \
    ${dir}/\${cluster_base}.nodes.tsv \
    ${dir}/\${cluster_base}.edges.tsv \
    $zfa_annotation_file
    """
}

process GBA_SUMMARY {
    label 'retry'
    publishDir "results", pattern: "plots/*.pdf"

    input:
    path("expts.txt")
    path("*")
    path("*")

    output:
    path("plots/*.pdf")

    script:
    """
    module load R/${params.r_version}
    gba-analysis.R --samples_file $params.expts
    """
}

process ENRICHMENT {
    label 'retry'
    publishDir "results", pattern: "*/GO/*"

    input:
    tuple val(dir), path(cluster_file), path(nodes_file)
    path(go_annotation_file)

    output:
    tuple val(dir), path(nodes_file), path("*/GO/*")

    script:
    """
    mkdir -p ${dir}/GO
    # run go enrichment script
    cluster_base=\$( basename $cluster_file )
    
    module load Python/$params.python_version
    create_files_for_topgo.py \
    --min_cluster_size $params.go_min_cluster_size \
    $nodes_file \$cluster_base

    if [[ \$( find ./ -type f -name "\${cluster_base}.cluster-*" | wc -l ) -le 1 ]] ; then
        echo "No clusters to test!"
        echo "No clusters to test!" > ${dir}/GO/done.txt
    else
        module load R/${params.r_version}
        module load topgo-wrapper/$params.ensembl_versionGO
        for file in \${cluster_base}.cluster-*
        do
            cluster_out=\$( basename -s '.tsv' \$file )
            echo \$cluster_out
            run_topgo.pl --input_file \${cluster_base}.all.tsv \
            --genes_of_interest_file \$file \
            --dir ${dir}/GO/\$cluster_out \
            --gene_field 3 \
            --name_field 4 \
            --description_field 5 \
            --go_terms_file $go_annotation_file
        done
    fi
    """
}

workflow {
    // Subset counts file to expts
    subset_output_ch = SUBSET(params.expts, params.all_counts)
    expts_file = subset_output_ch
        .flatMap { it[0] }
    sample_files = subset_output_ch
        .flatMap { it[1] }
        .map { [ it.parent.baseName, it ]}
    count_files = subset_output_ch
        .flatMap { it[2] }
        .map { [ it.parent.baseName, it ]}
    files_by_expt = sample_files.join(count_files)

    if ( params.debug ) {
        subset_output_ch.view { x -> "Subset output: $x" }
        expts_file.view { x -> "Expts file: $x"}
        sample_files.view { x -> "Expt name + sample file: $x"}
        count_files.view { x -> "Expt name + count file: $x"}
        files_by_expt.view { x -> "Expt name, sample and count files: $x" }
    }

    // Create a network for each expt
    orig_ch = CREATE_BASE_NETWORK(files_by_expt)
    expt_tab_ch = orig_ch.map { [it[0], it[2]]}
    cor_hist_ch = orig_ch.map { it[6] }
        .collect()
    mci_ch = orig_ch
        .map { [it[0], it[4]] }
    
    if ( params.debug ) {
        orig_ch.view { x -> "Base network output files: $x" }
        cor_hist_ch.view { x -> "Correlation histogram files: $x"}
        mci_ch.view { x -> "Expt name + mcx file: $x" }
    }

    // Test a range of filtering parameters
    stats_ch = TEST_PARAMETERS(mci_ch)
        .map { [it[1], it[2]] }
        .collect()
    if ( params.debug ) {
        stats_ch.view { x -> "Stats files: $x" }
    }

    // Filter networks by Correlation threshold or knn or both
    if (params.threshold) {
        threshold_params_ch = channel.value(params.threshold)
        filter_cor_ch = FILTER_COR(mci_ch, threshold_params_ch)
        if ( params.debug ) {
            filter_cor_ch.view { x -> "Filtered network files: $x" }
        }
    }
    if (params.knn) {
        knn_params_ch = channel.value(params.knn)
        filter_knn_ch = FILTER_KNN(mci_ch, knn_params_ch)
        if ( params.debug ) {
            filter_knn_ch.view { x -> "Filtered network files: $x" }
        }
    }
    if (params.threshold && params.knn) {
        filtered_ch = filter_cor_ch.concat(filter_knn_ch)
    } else if (params.threshold) {
        filtered_ch = filter_cor_ch
    } else {
        filtered_ch = filter_knn_ch
    }

    if ( params.debug ) {
        filtered_ch.view { x -> "MCI files to cluster: $x" }
    }

    filtered_stats_ch = filtered_ch.map { it[2] }
        .collect()
    
    // Collect up filtering stats and plot some graphs
    FILTER_STATS(expts_file, cor_hist_ch, stats_ch, filtered_stats_ch)
    if ( params.debug ) {
        filtered_stats_ch.view { x -> "Stats files from Filter processes: $x" }
    }

    if (params.clustering) {
        // Get annotation if necessary
        GET_ANNO_GET_GO_ANNO(
            params.species,
            [ anno_file: file(params.annotation_file), 
            version: params.ensembl_version,
            anno_script_url: params.get_anno_script_url,
            anno_bash_url: params.get_anno_bash_url ],
            [ go_anno_file: file(params.go_annotation_file),
            go_version: params.ensembl_versionGO,
            go_bash_url: params.get_go_anno_bash_url ]
        )
        // Cluster filtered networks
        infl_values_ch = channel.value(params.inflation_params)
        cluster_ch = CLUSTER(filtered_ch, infl_values_ch)

        // Run Guilt-by-Association on clustered networks
        tab_ch = expt_tab_ch.cross(cluster_ch)
            .map { [it[0][0], it[0][1], it[1][1], it[1][2] ] }
        annotation_ch = GET_ANNO_GET_GO_ANNO.out.anno_file
        go_file_ch = GET_ANNO_GET_GO_ANNO.out.go_anno_file
        zfa_annotation_ch = channel.value(params.zfa_file)
        gba_ch = GBA(tab_ch, annotation_ch, go_file_ch,
                       zfa_annotation_ch)

        if ( params.debug ) {
            annotation_ch.view { x -> "Annotation file: $x" }
            go_file_ch.view { x -> "GO annotation file: $x" }
            cluster_ch.view { x -> "Clustered MCI file: $x" }
            tab_ch.view { x -> "Tab file with clustered MCI file: $x" }
            gba_ch.view { x -> "Graph output files: $x" }
        }

        // Run GO enrichment on the clusters from the networks
        nodes_ch = gba_ch
            .map( { [it[0], it[1], it[3]] } )
        enrich_ch = ENRICHMENT(nodes_ch, go_file_ch)

        // Collect GBA stats and plot some graphs
        auc_stats_ch = gba_ch
            .map( { it[5] } )
            .collect()
        cluster_sizes_ch = cluster_ch
            .map { it[4][3] }
            .collect()
        gba_summary_ch = GBA_SUMMARY(expts_file, auc_stats_ch, cluster_sizes_ch)

        if ( params.debug ) {
            auc_stats_ch.view { x -> "AUC files: $x" }
            cluster_sizes_ch.view { x -> "Cluster sizes files: $x" }
            gba_summary_ch.view { x -> "GBA summary plot files: $x" }
        }
    }
}
