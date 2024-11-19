#!/usr/bin/env nextflow

log.info """\
  NETWORK CONSTRUCTION PIPELINE
  -----------------------------

  Testing: ${params.testing}
  Clustering: ${params.clustering}
  Expt to sample file: ${params.expts}
  All counts file: ${params.all_counts}
  Threshold params: ${params.threshold}
  KNN params: ${params.knn}
  Inflation Params: ${params.inflationParams}
"""

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
    """
    $params.QsubDir/subset-by-expt.sh -s ${params.ScriptDir} \
    ${expt_file} ${all_counts_file}
"""
}

process CREATE_BASE_NETWORK {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/*all-tpm*"

    input: 
    tuple val(expt_dir), path(sample_file), path(count_file)

    output:
    tuple val(expt_dir),
        path("$expt_dir/all-tpm.tsv"), path("$expt_dir/all-tpm.tab"),
        path("$expt_dir/all-tpm-orig.mcx"), path("$expt_dir/all-tpm-t20.mcx"),
        path("$expt_dir/all-tpm-orig.mat.csv.gz"),
        path("$expt_dir/$expt_dir-all-tpm-orig.cor-hist.txt")

    script:
    """
    mkdir -p $expt_dir
    awk -F"\\t" '{if(NR > 1){ print \$2 "\\t" \$3 }}' \
     $sample_file > $expt_dir/samples.txt

    module load R/$params.RVersion
    Rscript $params.ScriptDir/counts-to-fpkm-tpm.R \
    --transcripts_file $params.transcriptFile \
    --output_base $expt_dir/all --output_format tsv \
    --tpm $expt_dir/samples.txt $count_file

    module load MCL/$params.mclVersion

    # make network with all edges in
    mcxarray -data $expt_dir/all-tpm.tsv -co 0 \
    $params.skipRows $params.skipCols \
    $params.corMeasure $params.labels \
    --write-binary -o $expt_dir/all-tpm-orig.mcx \
    -write-tab $expt_dir/all-tpm.tab

    mcxdump -imx $expt_dir/all-tpm-orig.mcx \
    -tab $expt_dir/all-tpm.tab --dump-table \
    -digits 3 -sep-field "," -sep-lead "," \
    -o $expt_dir/all-tpm-orig.mat.csv
    gzip $expt_dir/all-tpm-orig.mat.csv

    module load Python/$params.PythonVersion
    python ${params.ScriptDir}/cor-hist.py \
    $expt_dir/all-tpm-orig.mat.csv.gz $expt_dir/$expt_dir-all-tpm-orig.cor-hist.txt

    # Also make one filtered with abs(), gt > 0.2
    mcx alter -imx $expt_dir/all-tpm-orig.mcx \
    -tf "abs(), gt(0.2)" \
    --write-binary -o $expt_dir/all-tpm-t20.mcx
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
    module load MCL/$params.mclVersion
    
    mkdir -p $dir

    # vary correlation
    mcx query -imx ${mci_file} --vary-correlation \
    --output-table > $dir/$dir-all-tpm.cor-stats.tsv

    # test varying k-nearest neighbours
    mcx query -imx ${mci_file} -vary-knn $params.knnTestParams \
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
    module load MCL/$params.mclVersion

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
    module load MCL/$params.mclVersion

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
    module load R/$params.RVersion
    Rscript $params.ScriptDir/edge-filtering-analysis.R \
    --samples_file $params.expts
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
    mkdir ${dir}

    module load MCL/${params.mclVersion}
    mcl $mci_file -I $inflation -o $dir/${mci_file}.I${inflationSuffix}
    
    mci_base=\$( basename $mci_file .mci )
    clm info $mci_file $dir/${mci_file}.I${inflationSuffix} >> $dir/\${mci_base}.info.txt

    clm info --node-all-measures --node-self-measures $mci_file \
    $dir/${mci_file}.I${inflationSuffix} > $dir/${mci_file}.I${inflationSuffix}.stats.tsv

    module load Python/$params.PythonVersion
    python ${params.ScriptDir}/summarise_clustering.py \
    --expt_name $dir $dir/${mci_file}.I${inflationSuffix}
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
        path("*/all-tpm*.auc.tsv"), path("*/all-tpm*.gene-scores.tsv"),
        path("*/all-tpm*.GBA-plots.pdf")

    script:
    matches = (mci_file =~ /all-tpm-(t?)(\d*)-?(k?)(\d*).mci$/)
    t_num = get_threshold(matches)
    println t_num
    """
    mkdir -p ${dir}
    cluster_base=\$( basename $cluster_file )
    
    # run convert_mcl script
    module load Python/$params.PythonVersion
    python ${params.ScriptDir}/convert_mcl.py \
    --min_cluster_size 4 --graph_id \${cluster_base} \
    --graphml_file ${dir}/\${cluster_base}.graphml \
    --nodes_file ${dir}/\${cluster_base}.nodes.tsv \
    --edges_file ${dir}/\${cluster_base}.edges.tsv \
    --edge_offset ${t_num} \
    $mci_file $cluster_file $tab_file $annotation_file

    module load R/$params.RVersion
    Rscript ${params.ScriptDir}/run-GBA-network.R \
    --auc_file ${dir}/\${cluster_base}.go.auc.tsv \
    --scores_file ${dir}/\${cluster_base}.go.gene-scores.tsv \
    --plots_file ${dir}/\${cluster_base}.go.GBA-plots.pdf \
    --min.term.size $params.minTermSize --max.term.size $params.maxTermSize \
    ${dir}/\${cluster_base}.nodes.tsv \
    ${dir}/\${cluster_base}.edges.tsv \
    $go_annotation_file

    Rscript ${params.ScriptDir}/run-GBA-network.R \
    --auc_file ${dir}/\${cluster_base}.zfa.auc.tsv \
    --scores_file ${dir}/\${cluster_base}.zfa.gene-scores.tsv \
    --plots_file ${dir}/\${cluster_base}.zfa.GBA-plots.pdf \
    --min.term.size $params.minTermSize --max.term.size $params.maxTermSize \
    ${dir}/\${cluster_base}.nodes.tsv \
    ${dir}/\${cluster_base}.edges.tsv \
    $zfa_annotation_file
    """
}

// process ANNOTATION {
//     label 'retry'
//     publishDir "reference", pattern: "*"

//     input:
//     tuple val(go_url)

//     output:
//     tuple path(go_annotation)

//     script:
//     """
//     wget $params.GO_URL
//     """
// }

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
    
    module load Python/$params.PythonVersion
    python ${params.ScriptDir}/create_files_for_topgo.py \
    --min_cluster_size $params.goMinClusterSize \
    $nodes_file \$cluster_base

    if [[ \$( find ./ -type f -name "\${cluster_base}.cluster-*" | wc -l ) -eq 0 ]] ; then
        echo "No clusters to test!"
        echo "No clusters to test!" > ${dir}/GO/done.txt
    else
        module load R/$params.RVersion
        module load topgo-wrapper/$params.EnsemblVersion
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

    subset_output_ch = SUBSET(params.expts, params.all_counts)
        .view()
    sample_files = subset_output_ch
        .flatMap { it[1] }
        .map { [ it.parent.baseName, it ]}
        .view()
    count_files = subset_output_ch
        .flatMap { it[2] }
        .map { [ it.parent.baseName, it ]}
        .view()
    files_by_expt = sample_files.join(count_files).view()

    orig_ch = CREATE_BASE_NETWORK(files_by_expt).view()
    expt_tab_ch = orig_ch.map { [it[0], it[2]]}

    tpm_ch = orig_ch
        .map { [it[0], it[1]] }
        .view()

    test_params_ch = TEST_PARAMETERS(tpm_ch)
        .map { [it[0], it[1]] }
        .view()

    if (params.threshold) {
        threshold_params_ch = channel.value(params.threshold)
        filter_cor_ch = FILTER_COR(test_params_ch, threshold_params_ch)
            .view()
    }

    if (params.knn) {
        knn_params_ch = channel.value(params.knn)
        filter_knn_ch = FILTER_KNN(test_params_ch, knn_params_ch)
            .view()
    }

    if (params.clustering) {

        if (params.threshold && params.knn) {
            filtered_ch = filter_cor_ch.concat(filter_knn_ch)
                .view()
        } else if (params.threshold) {
            filtered_ch = filter_cor_ch
        } else {
            filtered_ch = filter_knn_ch
        }
        infl_values_ch = channel.value(params.inflationParams)
        cluster_ch = CLUSTER(filtered_ch, infl_values_ch)
            .view()

        mci_ch = expt_tab_ch.cross(cluster_ch)
            .map { [it[0][0], it[0][1], it[1][1], it[1][2] ] }
            .view()
        
        annotation_ch = channel.value(params.AnnotationFile)
        go_annotation_ch = channel.value(params.GOFile)
        zfa_annotation_ch = channel.value(params.ZFAFile)
        gba_ch = GBA(mci_ch, annotation_ch, go_annotation_ch,
                       zfa_annotation_ch)
            .view()

        nodes_ch = gba_ch
            .map( { [it[0], it[1], it[3]] } )
        enrich_ch = ENRICHMENT(nodes_ch, go_annotation_ch)
    }
}
