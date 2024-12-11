#!/usr/bin/env nextflow

log.info """\
  NETWORK CONSTRUCTION PIPELINE
  -----------------------------

  Clustering: ${params.clustering}
  Debug: ${params.debug}
  Expt to sample file: ${params.samples}
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
process SUBSET_COUNTS {
    label 'process_high'
    publishDir "results", pattern: "expts.txt"

    input: 
    val sample_file
    val all_counts_file

    output:
    path("expts.txt"),               emit: expts_file
    path("*/samples.tsv"),           emit: sample_files // tuple of sample file names
    path("*/counts-by-gene.tsv"),    emit: count_files  // tuple of count file names

    script:
    debug = params.debug ? "-d" : ""
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

    input: 
    tuple val(expt), path(sample_file), path(count_file)

    output:
    tuple val(expt), path("${expt}-all-tpm.tsv"),       emit: tpms_file
    tuple val(expt), path("${expt}-all-tpm.tab"),       emit: tab_file // mapping of node ids to Ensembl gene ids
    tuple val(expt), path("${expt}-all-tpm-orig.mcx"),  emit: base_network
    tuple val(expt), path("${expt}-all-tpm-t20.mcx"),   emit: filtered_network
    path("${expt}-all-tpm-orig.cor-hist.txt"),          emit: cor_hist

    script:
    """
    awk -F"\\t" '{if(NR > 1){ print \$2 "\\t" \$3 }}' \
     ${sample_file} > ${expt}-samples.txt

    module load R/${params.r_version}
    counts-to-fpkm-tpm.R \
    --transcripts_file $params.transcript_file \
    --output_base ${expt}-all --output_format tsv \
    --tpm ${expt}-samples.txt $count_file

    module load MCL/$params.mcl_version

    # make network with all edges in
    mcxarray -data ${expt}-all-tpm.tsv -co 0 \
    $params.skip_rows $params.skip_cols \
    $params.cor_measure $params.labels \
    --write-binary -o ${expt}-all-tpm-orig.mcx \
    -write-tab ${expt}-all-tpm.tab

    # make histogram of correlation values
    mcx query -imx ${expt}-all-tpm-orig.mcx \
    -values-hist -1/1/20 \
    -o ${expt}-all-tpm-orig.cor-hist.txt

    # Also make one filtered with abs(), gt > 0.2
    mcx alter -imx ${expt}-all-tpm-orig.mcx \
    -tf "abs(), gt(0.2)" \
    --write-binary -o ${expt}-all-tpm-t20.mcx
    """

    stub:
    """
    mkdir -p ${expt}
    touch ${expt}-all-tpm.tsv ${expt}-all-tpm.tab \
    ${expt}-all-tpm-orig.mcx ${expt}-all-tpm-t20.mcx \
    ${expt}-all-tpm-orig.cor-hist.txt
    """
}

// Create a basic network with a correlation threshold of 0.2
// Test varying threshold and knn parameters
process TEST_PARAMETERS {
    label 'process_medium'

    input:
    tuple val(dir), path(mcx_file)

    output:
    tuple path("$dir-all-tpm.vary-cor-stats.tsv"),
        path("$dir-all-tpm.vary-knn-stats.tsv"),     emit: vary_threshold_stats

    script:
    """
    module load MCL/$params.mcl_version

    # vary correlation
    mcx query -imx ${mcx_file} --vary-correlation \
    --output-table > ${dir}-all-tpm.vary-cor-stats.tsv

    # test varying k-nearest neighbours
    mcx query -imx ${mcx_file} -vary-knn $params.knn_test_params \
    --output-table > ${dir}-all-tpm.vary-knn-stats.tsv
    """

    stub:
    """
    touch $dir-all-tpm.vary-cor-stats.tsv $dir-all-tpm.vary-knn-stats.tsv
    """
}

// Filter network by correlation threshold using MCL
process FILTER_COR {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mcl:14.137--0':
        'biocontainers/mcl:14.137--0' }"
    
    input:
    tuple val(dir), path(mcx_file)
    each threshold

    output:
    tuple val(dir), path("${dir}-all-tpm-t[0-9]*.mcx"),     emit: filtered_mci
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
    label 'process_low'

    input:
    tuple val(dir), path(mcx_file)
    each knn_threshold

    output:
    tuple val(dir), path("${dir}-all-tpm-t[0-9]*-k[0-9]*.mcx"), emit: filtered_mci
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
    label 'process_single'
    
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
        path("${dir}-all-tpm*.edges.tsv"),   emit: graph_files
    path("${dir}-all-tpm*.auc.tsv"),                                            emit: auc_files
    tuple path("${dir}-all-tpm*.gene-scores.tsv"), 
        path("${dir}-all-tpm*.GBA-plots.pdf"),                                  emit: gba_out

    script:
    matches = (mcx_file =~ /all-tpm-(t?)(\d*)-?(k?)(\d*).mcx$/)
    t_num = get_threshold(matches)
    if (params.debug){
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
    post-gba-analysis.R --samples_file ${sample_file}
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
    label 'process_single'

    input:
    tuple val(dir), path(cluster_file), path(graphml_file),
        path(nodes_file), path(edges_file)
    path(go_annotation_file)

    output:
    tuple val(dir), path(cluster_file), path("${cluster_file}.GO/*"), emit: go_output

    script:
    if (params.debug){
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
        echo "No clusters to test!" > ${cluster_file}.GO/done.txt
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

workflow {
    // Subset counts file to expts
    SUBSET_COUNTS(params.samples, params.all_counts)
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

    if ( params.debug ) {
        SUBSET_COUNTS.out.expts_file.view { x -> "Expts file: $x"}
        sample_files.view { x -> "Expt name + sample file: $x"}
        count_files.view { x -> "Expt name + count file: $x"}
        files_by_expt.view { x -> "Expt name, sample and count files: $x" }
    }

    // Create a network for each expt
    CREATE_BASE_NETWORK(files_by_expt)

    // Test a range of filtering parameters
    TEST_PARAMETERS(CREATE_BASE_NETWORK.out.filtered_network)

    // Filter networks by Correlation threshold or knn or both
    if (params.threshold) {
        threshold_params_ch = channel.value(params.threshold)
        FILTER_COR(CREATE_BASE_NETWORK.out.base_network, threshold_params_ch)
    }
    if (params.knn) {
        knn_params_ch = channel.value(params.knn)
        FILTER_KNN(CREATE_BASE_NETWORK.out.base_network, knn_params_ch)
    }

    // If both cor and knn have been run, concat the two channels together
    // Else the new channels are whichever one was run
    if (params.threshold && params.knn) {
        filtered_ch = FILTER_COR.out.filtered_mci.concat(FILTER_KNN.out.filtered_mci)
        filtered_stats_ch = FILTER_COR.out.node_stats
            .concat(FILTER_KNN.out.node_stats)
            .collect()
    } else if (params.threshold) {
        filtered_ch = FILTER_COR.out.filtered_mci
        filtered_stats_ch = FILTER_COR.out.node_stats.collect()
    } else {
        filtered_ch = FILTER_KNN.out.filtered_mci
        filtered_stats_ch = FILTER_KNN.out.node_stats.collect()
    }

    // Plot some graphs from the stats files
    // Collect together stats files
    cor_hist_ch = CREATE_BASE_NETWORK.out.cor_hist.collect()
    stats_ch = TEST_PARAMETERS.out.vary_threshold_stats.collect()
    FILTER_STATS(
        SUBSET_COUNTS.out.expts_file,
        params.samples,
        cor_hist_ch,
        stats_ch,
        filtered_stats_ch
    )
    if ( params.debug ) {
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
            [ anno_file: file(params.annotation_file), 
            version: params.ensembl_version,
            anno_script_url: params.get_anno_script_url,
            anno_bash_url: params.get_anno_bash_url ],
            [ go_anno_file: file(params.go_annotation_file),
            go_version: params.ensembl_versionGO,
            go_bash_url: params.get_go_anno_bash_url ]
        )
        if ( params.debug ) {
            GET_ANNO_GET_GO_ANNO.out.anno_file.view { x -> "Annotation file: $x" }
            GET_ANNO_GET_GO_ANNO.out.go_anno_file.view { x -> "GO annotation file: $x" }
        }

        // Cluster filtered networks
        infl_values_ch = channel.value(params.inflation_params)
        CLUSTER_NETWORK(
            filtered_ch,    // [ expt_name, filtered_mcx_file ]
            infl_values_ch  // inflation
        )
        if ( params.debug ) {
            CLUSTER_NETWORK.out.clustering.view { x -> "Clustered MCI file: $x" }
            CLUSTER_NETWORK.out.cluster_sizes.view { x -> "Cluster size files: $x" }
            CLUSTER_NETWORK.out.stats.view { x -> "Clustering stats files: $x" }
        }

        // join tab file to CLUSTER_NETWORK clustering output channel by expt name
        tab_ch = CREATE_BASE_NETWORK.out.tab_file
            .join(CLUSTER_NETWORK.out.clustering)
        // Run Guilt-by-Association on clustered networks
        RUN_GUILT_BY_ASSOCIATION(
            tab_ch,                                 // [ expt_name, tab_file, mcx_file, cluster_file ]
            GET_ANNO_GET_GO_ANNO.out.anno_file,     // [ gene_annotation_file ]
            GET_ANNO_GET_GO_ANNO.out.go_anno_file,  // [ GO_annotation_file ]
            channel.value(params.zfa_file)          // [ ZFA_annotation_file ]
        )

        RUN_POST_GBA_STATS(
            SUBSET_COUNTS.out.expts_file,
            params.samples,
            RUN_GUILT_BY_ASSOCIATION.out.auc_files.collect(),
            CLUSTER_NETWORK.out.cluster_sizes.collect()
        )

        // Run GO enrichment on the clusters from the networks
        // Uses the nodes file output from RUN_GUILT_BY_ASSOCIATION
        RUN_GO_ENRICHMENT(
            RUN_GUILT_BY_ASSOCIATION.out.graph_files,
            GET_ANNO_GET_GO_ANNO.out.go_anno_file
        )

        if ( params.debug ) {
            tab_ch.view { x -> "Tab file with clustered MCI file: $x" }
            RUN_GUILT_BY_ASSOCIATION.out.graph_files.view { x -> "Graph output files: $x" }
            RUN_GUILT_BY_ASSOCIATION.out.auc_files.collect().view { x -> "AUC output files: $x" }
            RUN_GUILT_BY_ASSOCIATION.out.gba_out.view { x -> "GBA output files: $x" }

            RUN_GO_ENRICHMENT.out.go_output.view { x -> "GO_ENRICHMENT output files: $x" }

            RUN_POST_GBA_STATS.out.plots.view { x -> "RUN_POST_GBA_STATS plots: $x" }
            RUN_POST_GBA_STATS.out.tsv.view { x -> "RUN_POST_GBA_STATS tsv files: $x" }
            RUN_POST_GBA_STATS.out.html.view { x -> "RUN_POST_GBA_STATS html files: $x" }
            // publish_ch = gba_ch
            //     .map { [ it[0], it[2], it[3], it[4], it[5], it[6], it[7] ] }
            //     .join(cluster_ch)
            //     .groupTuple()
            //     .view { x -> "Files to choose for publishing: $x" }
        }

    }
}
