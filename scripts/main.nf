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

    input: 
    val expt_file
    val all_counts_file

    output:
    path("${params.DirPrefix}-*")

    script:
    """
    $params.QsubDir/subset-by-expt.sh -s ${params.ScriptDir} \
-o ${params.DirPrefix} ${expt_file} ${all_counts_file}
"""
}

process CREATE_BASE_NETWORK {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"

    input: 
    path expt_dir

    output:
    tuple env(dir), path("*/all-tpm.tsv"), path("*/all-tpm.tab"), path("*/all-tpm-orig.mci")

    script:
    """
    dir=\$( basename $expt_dir )
    awk -F"\\t" '{if(NR > 1){ print \$2 "\\t" \$3 }}' \
     $expt_dir/samples.tsv > $expt_dir/samples.txt

    module load R/$params.RVersion
    Rscript $params.ScriptDir/counts-to-fpkm-tpm.R \
    --transcripts_file $params.transcriptFile \
    --output_base \$dir/all --output_format tsv \
    --tpm $expt_dir/samples.txt $expt_dir/counts-by-gene.tsv

    module load MCL/$params.mclVersion

    mcxarray -data \$dir/all-tpm.tsv -co 0 \
    $params.skipRows $params.skipCols \
    $params.corMeasure $params.labels \
    -o \$dir/all-tpm-orig.mci -write-tab \$dir/all-tpm.tab
    """
}

// Create a basic network with a correlation threshold of 0.2
// Test varying threshold and knn parameters
process TEST_PARAMETERS {
    label 'med_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(tpms_file)

    output:
    tuple val(dir), path("$dir/all-tpm-20.mci"), 
        path("$dir/all-tpm.cor-stats.tsv"), path("$dir/all-tpm.knn-stats.tsv")

    script:
    """
    module load MCL/$params.mclVersion
    
    mkdir -p $dir
    mcxarray -data $tpms_file -co 0.2 \
    $params.skipRows $params.skipCols -tf 'abs()' \
    $params.corMeasure $params.labels -o $dir/all-tpm-20.mci
    
    # vary correlation
    mcx query -imx $dir/all-tpm-20.mci --vary-correlation \
    --output-table > $dir/all-tpm.cor-stats.tsv

    # test varying k-nearest neighbours
    mcx query -imx $dir/all-tpm-20.mci -vary-knn $params.knnTestParams \
    --output-table > $dir/all-tpm.knn-stats.tsv
    """
}

// Threshold
process FILTER_COR {
    label 'retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(mci_file)
    each threshold

    output:
    tuple val(dir), path("$dir/*-[t][0-9]*.mci"),
        path("$dir/*stats.tsv")

    script:
    Integer suffix = threshold * 100
    outputBase = [dir, "/all-tpm"].join('')
    thresholdBase = [outputBase, "t$suffix"].join("-")
    """
    module load MCL/$params.mclVersion

    mkdir -p $dir
    mcx alter -imx ${mci_file} -tf \
    "gq(${threshold}), add(-${threshold})" \
    -o ${thresholdBase}.mci
    mcx query -imx ${thresholdBase}.mci > ${thresholdBase}.stats.tsv
    """
}

process FILTER_KNN {
    label 'retry'
    publishDir "results", pattern: "*/all-tpm*"

    input:
    tuple val(dir), path(mci_file)
    each knn_threshold

    output:
    tuple val(dir), path("$dir/*-[k][0-9]*.mci"),
        path("$dir/*stats.tsv")

    script:
    """
    module load MCL/$params.mclVersion

    mkdir -p $dir
    knnBase=${dir}/all-tpm-t20-k${knn_threshold}
    mcx alter -imx ${mci_file} -tf "add(-0.2), #knn($knn_threshold)" \
    -o \${knnBase}.mci
    mcx query -imx \${knnBase}.mci > \${knnBase}.stats.tsv
    """
}

// Cluster network
process CLUSTER {
    label 'retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(mci_file), path(stats_file)
    each inflation

    output:
    tuple val(dir), path(mci_file), path("*/*.mci.I[0-9][0-9]"),
        path("*/*.mci.I[0-9][0-9].stats.tsv"), path("*/*.mci.I[0-9][0-9].cl*")

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
    python ${params.ScriptDir}/summarise_clustering.py $dir/${mci_file}.I${inflationSuffix}
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

process ENRICHMENT {
    label 'retry'
    publishDir "results", pattern: "*/all-tpm*go*"

    input:
    tuple val(dir), path(cluster_file), path(nodes_file)

    output:
    tuple val(dir), path(nodes_file), path("*/all-tpm*go*")

    script:
    """
    mkdir -p ${dir}
    # run go enrichment script
    cluster_base=\$( basename $cluster_file )
    
    module load R/$params.RVersion
    Rscript ${params.ScriptDir}/gprofiler-on-network-clusters.R \
    --min_cluster_size ${params.goMinClusterSize} \
    --output_file ${dir}/\${cluster_base}.go-enrichments.tsv \
    --rds_file ${dir}/\${cluster_base}.go-enrichments.rds \
    $nodes_file
    """
}

workflow {

    subset_output_ch = SUBSET(params.expts, params.all_counts)
        .flatten()
        .view()

    orig_ch = CREATE_BASE_NETWORK(subset_output_ch).view()

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

        mci_ch = orig_ch.cross(cluster_ch)
            .map { [it[0][0], it[0][2], it[1][1], it[1][2] ] }
            .view()
        
        annotation_ch = channel.value(params.AnnotationFile)
        go_annotation_ch = channel.value(params.GOFile)
        zfa_annotation_ch = channel.value(params.ZFAFile)
        gba_ch = GBA(mci_ch, annotation_ch, go_annotation_ch,
                       zfa_annotation_ch)
            .view()

        nodes_ch = gba_ch
            .map( { [it[0], it[1], it[3]] } )
        enrich_ch = ENRICHMENT(nodes_ch)
            .view()
    }
}
