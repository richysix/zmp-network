#!/usr/bin/env nextflow

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

process LOG_INFO {
    label 'process_single'
    executor 'local'

    exec:
    log.info """\
    NETWORK CONSTRUCTION PIPELINE
    -----------------------------

    Clustering:             ${params.clustering}
    Debug:                  ${params.debug}
    Expt to sample file:    ${params.samples}
    All counts file:        ${params.all_counts}
    Reference dir:          ${params.ref_dir}
    Transcript File:        ${params.transcript_file}
    Annotation File:        ${params.annotation_file}
    Go Annotation File:     ${params.go_annotation_file}
    Threshold params:       ${params.threshold}
    KNN params:             ${params.knn}
    Inflation Params:       ${params.inflation_params}
    """

    stub:
    """
    echo "Log info"
    """
}

// process to subset the file containing all the samples
// to each separate experiment
// see README
process SUBSET_COUNTS {
    label 'process_medium'
    publishDir "results", pattern: "expts.txt"

    input: 
    path(sample_file)
    path(all_counts_file)

    output:
    path("expts.txt"),               emit: expts_file
    path("*/samples.tsv"),           emit: sample_files // tuple of sample file names
    path("*/counts-by-gene.tsv"),    emit: count_files  // tuple of count file names

    script:
    debug = params.debug > 0 ? "-d" : ""
    """
    subset-by-expt.sh ${debug} ${sample_file} ${all_counts_file}
    """

    stub:
    """
    touch expts.txt
    for dir in test-1 test-2
    do
      mkdir \$dir
      touch \$dir/samples.tsv
      touch \$dir/counts-by-gene.tsv
    done
    """
}

// Create correlation network containing all edges
// Also create one, removing edges abs() < 0.2 for testing cut-off parameters
// in TEST_PARAMETERS (smaller to load into memory)
process CREATE_BASE_NETWORK {
    label 'process_medium'
    publishDir path: "${params.outdir}/${expt}", mode: params.publish_dir_mode,
        pattern: "${expt}-tpm*{-orig.mcx,-orig.cor-hist.txt,.tsv,.tab,}"


    input: 
    tuple val(expt), path(sample_file), path(count_file)
    path(transcript_file)

    output:
    tuple val(expt), path("${expt}-tpm.tsv"),                       emit: tpms_file
    tuple val(expt), path("${expt}-tpm.tab"),                       emit: tab_file // mapping of node ids to Ensembl gene ids
    tuple val(expt), path("${expt}-tpm-orig.mcx"),                  emit: base_network
    tuple val(expt), path("${expt}-tpm-filtered-by-zeros.tsv"),     emit: filtered_tpms_file
    tuple val(expt), path("${expt}-tpm-filtered.tab"),              emit: filtered_tab_file // mapping of node ids to Ensembl gene ids
    tuple val(expt), path("${expt}-tpm-filtered-orig.mcx"),         emit: filtered_network
    tuple val(expt), path("${expt}-tpm-filtered-t20.mcx"),          emit: filtered_t20_network
    path("${expt}-tpm-orig.cor-hist.txt"),                          emit: cor_hist
    path("${expt}-tpm-filtered-orig.cor-hist.txt"),                 emit: filtered_cor_hist

    script:
    """
    awk -F"\\t" '{if(NR > 1){ print \$2 "\\t" \$3 }}' \
        ${sample_file} > ${expt}-samples.txt

    module load R/${params.r_version}
    counts-to-fpkm-tpm.R \
    --transcripts_file ${transcript_file} \
    --output_base ${expt} --output_format tsv \
    --tpm --filter_threshold ${params.high_cor_filter_threshold} \
    ${expt}-samples.txt ${count_file}

    module load MCL/$params.mcl_version

    # make network with all genes and edges in
    mcxarray -data ${expt}-tpm.tsv -co 0 \
    $params.skip_rows $params.skip_cols \
    $params.cor_measure $params.labels \
    --write-binary -o ${expt}-tpm-orig.mcx \
    -write-tab ${expt}-tpm.tab

    # make histogram of correlation values
    mcx query -imx ${expt}-tpm-orig.mcx \
    -values-hist -1/1/20 \
    -o ${expt}-tpm-orig.cor-hist.txt

    # make network with filtered gene set, all edges
    mcxarray -data ${expt}-tpm-filtered-by-zeros.tsv -co 0 \
    $params.skip_rows $params.skip_cols \
    $params.cor_measure $params.labels \
    --write-binary -o ${expt}-tpm-filtered-orig.mcx \
    -write-tab ${expt}-tpm-filtered.tab

    # make histogram of correlation values
    mcx query -imx ${expt}-tpm-filtered-orig.mcx \
    -values-hist -1/1/20 \
    -o ${expt}-tpm-filtered-orig.cor-hist.txt

    # Also make one filtered with abs(), gt > 0.2
    mcx alter -imx ${expt}-tpm-filtered-orig.mcx \
    -tf "abs(), gt(0.2)" \
    --write-binary -o ${expt}-tpm-filtered-t20.mcx
    """

    stub:
    """
    touch ${expt}-samples.txt ${expt}-tpm.tsv ${expt}-tpm.tab \
    ${expt}-tpm-orig.mcx ${expt}-tpm-orig.cor-hist.txt \
    ${expt}-tpm-filtered-by-zeros.tsv ${expt}-tpm-filtered.tab \
    ${expt}-tpm-filtered-orig.mcx ${expt}-tpm-filtered-t20.mcx \
    ${expt}-tpm-filtered-orig.cor-hist.txt
    """
}

// Test varying threshold and knn parameters
process TEST_PARAMETERS {
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mcl:14.137--0':
        'biocontainers/mcl:14.137--0' }"

    input:
    tuple val(dir), path(mcx_file)

    output:
    tuple path("${mcx_base}.vary-cor-stats.tsv"),
        path("${mcx_base}.vary-knn-stats.tsv"),     emit: vary_threshold_stats

    script:
    mcx_base=mcx_file.baseName
    """
    # vary correlation
    mcx query -imx ${mcx_file} --vary-correlation \
    --output-table > ${mcx_base}.vary-cor-stats.tsv

    # test varying k-nearest neighbours
    mcx query -imx ${mcx_file} -vary-knn $params.knn_test_params \
    --output-table > ${mcx_base}.vary-knn-stats.tsv
    """

    stub:
    """
    touch ${mcx_base}.vary-cor-stats.tsv ${mcx_base}.vary-knn-stats.tsv
    """
}

// Filter network by correlation threshold using MCL
process FILTER_COR {
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mcl:14.137--0':
        'biocontainers/mcl:14.137--0' }"
    
    input:
    tuple val(dir), path(mcx_file)
    each threshold

    output:
    tuple val(dir), path("${dir}-all-tpm-t[0-9]*.mcx"),     emit: filtered_mcx
    path("${dir}-all-tpm-t[0-9]*.stats.tsv"),               emit: node_stats

    script:
    Integer suffix = threshold * 100
    """
    mcx alter -imx ${mcx_file} \
    -tf "gt(${threshold}), add(-${threshold})" \
    --write-binary -o ${dir}-all-tpm-t${suffix}.mcx
    mcx alter -imx ${dir}-all-tpm-t${suffix}.mcx \
    -tf "add(${threshold})" | mcx query -imx - > \
    ${dir}-all-tpm-t${suffix}.stats.tsv
    """

    stub:
    Integer suffix = threshold * 100
    def mcx = "${dir}-all-tpm-t${suffix}.mcx"
    def stats = "${dir}-all-tpm-t${suffix}.stats.tsv"
    """
    touch ${mcx} ${stats}
    """
}

// Filter network by K-nearest-neighbours using MCL
process FILTER_KNN {
    label 'process_medium'

    input:
    tuple val(dir), path(mcx_file)
    each knn_threshold

    output:
    tuple val(dir), path("${dir}-all-tpm-t[0-9]*-k[0-9]*.mcx"), emit: filtered_mcx
    path("${dir}-all-tpm-t[0-9]*-k[0-9]*.stats.tsv"),           emit: node_stats

    script:
    threshold=0.2
    thresholdSuffix=20
    """
    module load MCL/$params.mcl_version

    knnBase="${dir}-all-tpm-t${thresholdSuffix}-k${knn_threshold}"
    mcx alter -imx ${mcx_file} \
    -tf "abs(), gt(${threshold}), add(-${threshold}), #knn($knn_threshold)" \
    --write-binary -o \${knnBase}.mcx
    # stats
    mcx alter -imx \${knnBase}.mcx -tf "add(${threshold})" | \
    mcx query -imx - > \${knnBase}.stats.tsv
    """

    stub:
    def mci = "${dir}-all-tpm-t20-k${knn_threshold}.mcx"
    def stats = "${dir}-all-tpm-t20-k${knn_threshold}.stats.tsv"
    """
    touch ${mci} ${stats}
    """
}

// Produce some plots after filtering
process FILTER_STATS {
    label 'process_single'
    publishDir "results", pattern: "plots/*.pdf"

    input:
    path("expts.txt")   // list of expt names
    path(sample_file)   // samples file
    path("*")           // correlation histogram files
    path("*")           // vary threshold stats files
    path("*")           // node stats after filtering

    output:
    path("plots/*.pdf"), emit: plots

    script:
    """
    module load R/${params.r_version}
    edge-filtering-analysis.R --samples_file $sample_file
    """

    stub:
    """
    mkdir -p plots/
    touch plots/all-cor-dist.pdf plots/cor-stats-node-degrees-close-up.pdf \
    plots/knn-stats-all-components.pdf
    for method in cor knn 
    do
      touch plots/\$method-stats-singletons.pdf plots/\$method-stats-nodes.pdf \
      plots/\$method-stats-nodes-singletons.pdf \
      plots/\$method-stats-node-degrees.pdf \
      plots/\$method-stats-node-degree-distribution.pdf \
      plots/\$method-stats-node-degree-distribution-close-up.pdf
    done
    """
}

// Cluster network
process CLUSTER_NETWORK {
    label 'process_low'
    
    input:
    tuple val(dir), path(mcx_file)
    each inflation

    output:
    tuple val(dir), path(mcx_file), path("${mcx_file}.I[0-9][0-9]"),     emit: clustering
    path("${mcx_file}.I[0-9][0-9].cl-sizes.tsv"),                        emit: cluster_sizes
    tuple path("${mcx_file}.I[0-9][0-9].stats.tsv"),
        path("${mcx_file}.I[0-9][0-9].info.txt"),                        emit: stats

    script:
    Integer inflationSuffix = inflation * 10
    def file_base = "${mcx_file}.I${inflationSuffix}"
    """
    module load MCL/${params.mcl_version}
    mcl ${mcx_file} -I ${inflation} -o ${file_base}

    clm info ${mcx_file} ${file_base} >> ${file_base}.info.txt

    clm info --node-all-measures --node-self-measures ${mcx_file} \
    ${file_base} > ${file_base}.stats.tsv

    module load Python/${params.python_version}
    summarise_clustering.py ${file_base}
    """

    stub:
    Integer inflationSuffix = inflation * 10
    def mcx_cluster_file = "${mcx_file}.I${inflationSuffix}"
    """
    touch ${mcx_cluster_file} ${mcx_cluster_file}.stats.tsv \
    ${mcx_cluster_file}.cl-sizes.tsv ${mcx_cluster_file}.info.txt
    """
}

// Run Guilt by association analysis on the clustered network
// with both GO and ZFA annotation
process RUN_GUILT_BY_ASSOCIATION {
    label 'process_high_memory'

    input:
    tuple val(dir), path(tab_file), path(mcx_file), path(cluster_file)
    path(annotation_file)
    path(go_annotation_file)
    path(zfa_annotation_file)

    output:
    tuple path(cluster_file),
        path("${dir}-all-tpm*.graphml"),
        path("${dir}-all-tpm*.nodes.tsv"),
        path("${dir}-all-tpm*.edges.tsv"),              emit: graph_files
    path("${dir}-all-tpm*.auc.tsv"),                    emit: auc_files
    tuple path("${dir}-all-tpm*.gene-scores.tsv"), 
        path("${dir}-all-tpm*.GBA-plots.pdf"),          emit: gba_out

    script:
    matches = (mcx_file =~ /all-tpm-(t?)(\d*)-?(k?)(\d*).mcx$/)
    t_num = get_threshold(matches)
    if (params.debug > 1 ){
        println("GBA: Threshold value = " + t_num)
    }
    """
    # convert from binary
    module load MCL/${params.mcl_version}
    mcx_base=\$( basename ${mcx_file} .mcx)
    mcx convert ${mcx_file} \${mcx_base}.mci

    # run convert_mcl script
    module load Python/${params.python_version}
    convert_mcl.py \
    --min_cluster_size 4 --graph_id ${cluster_file} \
    --graphml_file ${cluster_file}.graphml \
    --nodes_file ${cluster_file}.nodes.tsv \
    --edges_file ${cluster_file}.edges.tsv \
    --edge_offset ${t_num} \
    \${mcx_base}.mci ${cluster_file} ${tab_file} ${annotation_file}

    module load R/${params.r_version}
    run-GBA-network.R \
    --auc_file ${cluster_file}.go.auc.tsv \
    --scores_file ${cluster_file}.go.gene-scores.tsv \
    --plots_file ${cluster_file}.go.GBA-plots.pdf \
    --min.term.size ${params.min_term_size} --max.term.size ${params.max_term_size} \
    ${cluster_file}.nodes.tsv \
    ${cluster_file}.edges.tsv \
    ${go_annotation_file}

    run-GBA-network.R \
    --auc_file ${cluster_file}.zfa.auc.tsv \
    --scores_file ${cluster_file}.zfa.gene-scores.tsv \
    --plots_file ${cluster_file}.zfa.GBA-plots.pdf \
    --min.term.size ${params.min_term_size} --max.term.size ${params.max_term_size} \
    ${cluster_file}.nodes.tsv \
    ${cluster_file}.edges.tsv \
    ${zfa_annotation_file}
    """

    stub:
    """
    for name in edges.tsv go.GBA-plots.pdf go.gene-scores.tsv graphml \
    nodes.tsv zfa.GBA-plots.pdf zfa.gene-scores.tsv
    do
        touch ${cluster_file}.\$name
    done
    touch ${cluster_file}.auc.tsv
    """
}

process RUN_POST_GBA_STATS {
    label 'process_single'
    publishDir "results", pattern: "plots/*.pdf"
    publishDir "results", pattern: "*.{tsv,html}"

    input:
    path("expts.txt")   // list of expt names
    path(sample_file)   // samples file
    path("*")           // GBA AUC output files
    path("*")           // CLUSTER cluster size files

    output:
    path("plots/*.pdf"),                    emit: plots
    tuple path("top_ecdf_diff.tsv"),
        path("cor-cluster_ecdf_diff.tsv"),
        path("knn-cluster_ecdf_diff.tsv"),  emit: tsv
    path("top_ecdf_diff_table-*.html"),     emit: html

    script:
    """
    module load R/${params.r_version}
    post-gba-analysis.R --samples_file ${sample_file} \
    --min_cluster_size 100 \
    --max_cluster_size 1000
    """

    stub:
    """
    mkdir -p plots
    touch \
      top_ecdf_diff.tsv \
      cor-cluster_ecdf_diff.tsv \
      knn-cluster_ecdf_diff.tsv \
      top_ecdf_diff_table-100-1000.html \
      plots/auc-diff.pdf \
      plots/cor-clustering-summary-plots.pdf \
      plots/knn-clustering-summary-plots.pdf
    """
}

process RUN_GO_ENRICHMENT {
    label 'process_long'

    input:
    tuple path(cluster_file), path(graphml_file),
        path(nodes_file), path(edges_file)
    path(go_annotation_file)

    output:
    tuple path(cluster_file), path("${cluster_file}.GO"), emit: go_output

    script:
    if (params.debug > 1 ){
        println("RUN_GO_ENRICHMENT: Name of cluster file is ${cluster_file}")
    }
    """
    mkdir -p ${cluster_file}.GO
    
    # run script to create input files
    module load Python/${params.python_version}
    create_files_for_topgo.py \
    --min_cluster_size ${params.go_min_cluster_size} \
    ${nodes_file} ${cluster_file}

    # run go enrichment script
    if [[ \$( find ./ -type f -name "${cluster_file}.cluster-*" | wc -l ) -le 1 ]] ; then
        echo "No clusters to test!"
        echo "No clusters to test!" > ${cluster_file}.GO/${cluster_file}-done.txt
    else
        module load R/${params.r_version}
        module load topgo-wrapper/${params.ensembl_versionGO}
        for file in ${cluster_file}.cluster-*
        do
            cluster_out=\$( basename -s '.tsv' \$file )
            echo \$cluster_out
            run_topgo.pl --input_file ${cluster_file}.all.tsv \
            --genes_of_interest_file \$file \
            --dir ${cluster_file}.GO/\$cluster_out \
            --gene_field 3 \
            --name_field 4 \
            --description_field 5 \
            --go_terms_file $go_annotation_file
        done
    fi
    """

    stub:
    """
    mkdir -p ${cluster_file}.GO
    for cluster in \$( seq 5 )
    do
      go_dir=${cluster_file}.GO/${cluster_file}.cluster-\$cluster
      mkdir \$go_dir
      for domain in BP MF CC
      do
        touch \$go_dir/\$domain.all.genes.tsv \$go_dir/\$domain.all.tsv \
        \$go_dir/\$domain.pdf \$go_dir/\$domain.sig.genes.tsv \
        \$go_dir/\$domain.sig.tsv
      done
    done
    """
}

process PUBLISH_NETWORKS {
    label 'process_single'
    publishDir path: params.outdir, mode: params.publish_dir_mode,
        pattern: "${expt}/*"

    input:
    tuple val(expt), val(method) // name of experiment
    path("*") // RUN_POST_GBA_STATS.out.tsv: ecdf tsv files
    // path("*") // base_network_ch: tpm, tab and base network files
    path("*") // filtered_network_ch: t*-k*-mcx
    path("*") // stats_files_ch: network stats files
    path("*") // graph_files_with_go_ch: graphs and GO output
    path("*") // auc_files_ch: AUC files
    path("*") // RUN_GUILT_BY_ASSOCIATION.out.gba_out.flatten().collect(): other GBA output files

    output:
    path("${expt}/*"), emit: pub_files

    script:
    """
    module load Python/${params.python_version}
    stage_output_files.py --expt ${expt} --method ${method} top_ecdf_diff.tsv
    """

    stub:
    """
    for dir in test-1 test-2
    do
        mkdir -p \$dir
        base="all-tpm-t20-k3.mcx"
        mv \$dir-\$base \$dir/
        mv \$dir-\$base.I14 \$dir/
        for suffix in cl-sizes.tsv stats.tsv \
        edges.tsv nodes.tsv graphml \
        go.GBA-plots.pdf zfa.GBA-plots.pdf \
        go.gene-scores.tsv zfa.gene-scores.tsv
        do
            mv \$dir-\$base.I14.\$suffix \$dir/
        done
        mv \$dir-\$base.I14.GO \$dir/
    done
    """
}

workflow {
    if ( params.debug > 0 ) {
        LOG_INFO()
    }
    // Subset counts file to expts
    SUBSET_COUNTS(
        file(params.samples),
        file(params.all_counts)
    )
    // extract expt name from path to use as join key
    // it.parent.baseName gets final directory name (e.g. test-1 from /work/78/810e9b6891dd6476b1474a90983952/test-1/samples.tsv)
    sample_files = SUBSET_COUNTS.out.sample_files
        .flatten()
        .map { [ it.parent.baseName, it ]}
    // extract expt name from path to use as join key
    count_files = SUBSET_COUNTS.out.count_files
        .flatten()
        .map { [ it.parent.baseName, it ]}
    // join sample files to count files by expt name
    files_by_expt = sample_files.join(count_files)

    if ( params.debug > 1 ) {
        SUBSET_COUNTS.out.expts_file.view { x -> "Expts file: $x"}
        sample_files.view { x -> "Expt name + sample file: $x"}
        count_files.view { x -> "Expt name + count file: $x"}
        files_by_expt.view { x -> "Expt name, sample and count files: $x" }
    }

    // Create a network for each expt
    CREATE_BASE_NETWORK(
        files_by_expt,
        file(params.transcript_file)
    )

    // Test a range of filtering parameters
    TEST_PARAMETERS(CREATE_BASE_NETWORK.out.filtered_t20_network)

    // Filter networks by correlation threshold or knn or both
    if (params.threshold) {
        threshold_params_ch = channel.value(params.threshold)
        FILTER_COR(CREATE_BASE_NETWORK.out.filtered_network, threshold_params_ch)
    }
    if (params.knn) {
        knn_params_ch = channel.value(params.knn)
        FILTER_KNN(CREATE_BASE_NETWORK.out.filtered_network, knn_params_ch)
    }

    // If both cor and knn have been run, concat the two channels together
    // Else the new channels are whichever one was run
    if (params.threshold && params.knn) {
        filtered_ch = FILTER_COR.out.filtered_mcx.concat(FILTER_KNN.out.filtered_mcx)
        filtered_stats_ch = FILTER_COR.out.node_stats
            .concat(FILTER_KNN.out.node_stats)
            .collect()
    } else if (params.threshold) {
        filtered_ch = FILTER_COR.out.filtered_mcx
        filtered_stats_ch = FILTER_COR.out.node_stats.collect()
    } else {
        filtered_ch = FILTER_KNN.out.filtered_mcx
        filtered_stats_ch = FILTER_KNN.out.node_stats.collect()
    }

    // Plot some graphs from the stats files
    // Collect together stats files
    cor_hist_ch = CREATE_BASE_NETWORK.out.cor_hist.collect()
    filtered_cor_hist_ch = CREATE_BASE_NETWORK.out.filtered_cor_hist.collect()
    stats_ch = TEST_PARAMETERS.out.vary_threshold_stats.collect()
    FILTER_STATS(
        SUBSET_COUNTS.out.expts_file,
        file(params.samples),
        cor_hist_ch,
        stats_ch,
        filtered_stats_ch
    )
    if ( params.debug > 1 ) {
        // Filtered network files for clustering
        filtered_ch.view { x -> "Filtered network files: $x" }
        // Files for FILTER_STATS
        cor_hist_ch.view { x -> "Cor histogram files: $x" }
        stats_ch.view { x -> "Vary threshold stats files: $x" }
        filtered_stats_ch.view { x -> "Node stats from filtering: $x" }
    }

    if (params.clustering) {
        // Get annotation if necessary
        GET_ANNO_GET_GO_ANNO(
            params.species,
            [ 
                anno_file: file(params.annotation_file), 
                version: params.ensembl_version
            // anno_script_url: params.get_anno_script_url,
            // anno_bash_url: params.get_anno_bash_url 
            ],
            [ 
                go_anno_file: file(params.go_annotation_file),
                go_version: params.ensembl_versionGO
            ]
        )
        if ( params.debug > 1 ) {
            GET_ANNO_GET_GO_ANNO.out.anno_file.view { x -> "Annotation file: $x" }
            GET_ANNO_GET_GO_ANNO.out.go_anno_file.view { x -> "GO annotation file: $x" }
        }

        // Cluster filtered networks
        infl_values_ch = channel.value(params.inflation_params)
        CLUSTER_NETWORK(
            filtered_ch,    // [ expt_name, filtered_mcx_file ]
            infl_values_ch  // inflation
        )
        if ( params.debug > 1 ) {
            CLUSTER_NETWORK.out.clustering.view { x -> "Clustered MCI file: $x" }
            CLUSTER_NETWORK.out.cluster_sizes.view { x -> "Cluster size files: $x" }
            CLUSTER_NETWORK.out.stats.view { x -> "Clustering stats files: $x" }
        }

        // join tab file to CLUSTER_NETWORK clustering output channel by expt name
        // Have to use the combine operator because the key (expt_name)
        // is not unique, because the networks are clustered with different inflation params
        tab_ch = CREATE_BASE_NETWORK.out.filtered_tab_file
            .combine(CLUSTER_NETWORK.out.clustering, by: 0)

        // Run Guilt-by-Association on clustered networks
        RUN_GUILT_BY_ASSOCIATION(
            tab_ch,                                 // [ expt_name, tab_file, mcx_file, cluster_file ]
            GET_ANNO_GET_GO_ANNO.out.anno_file,     // [ gene_annotation_file ]
            GET_ANNO_GET_GO_ANNO.out.go_anno_file,  // [ GO_annotation_file ]
            channel.value(file(params.zfa_file))    // [ ZFA_annotation_file ]
        )

        auc_files_ch = RUN_GUILT_BY_ASSOCIATION.out.auc_files.collect()
        RUN_POST_GBA_STATS(
            SUBSET_COUNTS.out.expts_file,
            file(params.samples),
            auc_files_ch,
            CLUSTER_NETWORK.out.cluster_sizes.collect()
        )

        // Run GO enrichment on the clusters from the networks
        // Uses the nodes file output from RUN_GUILT_BY_ASSOCIATION
        RUN_GO_ENRICHMENT(
            RUN_GUILT_BY_ASSOCIATION.out.graph_files,
            GET_ANNO_GET_GO_ANNO.out.go_anno_file
        )

        filtered_network_ch = filtered_ch.map { _expt, mcx -> mcx }
            .collect()
        // stats files from clustering
        stats_files_ch = CLUSTER_NETWORK.out.cluster_sizes
            .mix(CLUSTER_NETWORK.out.stats)
            .collect()

        // join graph files to GO enrichment out file by clustered network
        // Need to remove the path from the network files for the join to work
        graph_files_ch = RUN_GUILT_BY_ASSOCIATION.out.graph_files
            .map { network, graphml, node, edges -> [ network.name, network, graphml, node, edges ] }
        go_out_files_ch = RUN_GO_ENRICHMENT.out.go_output
            .map { network, go_out_files -> [ network.name, go_out_files ] }
        // join together with network name, the remove it so that everything is a path
        graph_files_with_go_ch = graph_files_ch.join(go_out_files_ch)
            .map { _name, mcx, graphml, node, edges, go_out_files -> [ mcx, graphml, node, edges, go_out_files ] }
            .flatten()
            .collect()

        expts_ch = sample_files.map { expt, _samples_file -> [ expt, "knn" ] }
        PUBLISH_NETWORKS(
            expts_ch,
            RUN_POST_GBA_STATS.out.tsv,
            // base_network_ch,
            filtered_network_ch,
            stats_files_ch,
            graph_files_with_go_ch,                         // [ mcx.I*, graphml, node, edges, GO_files]
            auc_files_ch,
            RUN_GUILT_BY_ASSOCIATION.out.gba_out.flatten().collect()
        )

        if ( params.debug > 1 ) {
            tab_ch.view { x -> "Tab file with clustered MCX file: $x" }
            RUN_GUILT_BY_ASSOCIATION.out.graph_files.view { x -> "Graph output files: $x" }
            RUN_GUILT_BY_ASSOCIATION.out.auc_files.collect().view { x -> "AUC output files: $x" }
            RUN_GUILT_BY_ASSOCIATION.out.gba_out.flatten().collect().view { x -> "GBA output files: $x" }

            RUN_GO_ENRICHMENT.out.go_output.view { x -> "GO_ENRICHMENT output files: $x" }

            RUN_POST_GBA_STATS.out.plots.view { x -> "RUN_POST_GBA_STATS plots: $x" }
            RUN_POST_GBA_STATS.out.tsv.view { x -> "RUN_POST_GBA_STATS tsv files: $x" }
            RUN_POST_GBA_STATS.out.html.view { x -> "RUN_POST_GBA_STATS html files: $x" }

            filtered_network_ch.view { x -> "filtered mcx: $x" }
            graph_files_with_go_ch.view { x -> "Graph files with GO enrichment: $x" }
        }

    }
}
