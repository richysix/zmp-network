#!/usr/bin/env nextflow

log.info """\
  NETWORK CONSTRUCTION PIPELINE
  -----------------------------

  Testing: ${params.testing}
  Expt to sample file: ${params.expts}
  All counts file: ${params.all_counts}
"""

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
    --transcripts_file $params.RefDir/$params.transcriptFile \
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
    label 'big_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(tpms_file)

    output:
    tuple val(dir), path("$dir/all-tpm-20.mci"), 
        path("$dir/all-tpm.cor-stats.tsv"), path("$dir/all-tpm.knn-stats.tsv")

    script:
    """
    module load MCL/$params.mclVersion
    
    mkdir $dir
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
process THRESHOLD {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(mci_file)

    output:
    tuple val(dir), path("$dir/*-[kt][0-9]*.mci"),
        path("$dir/*stats.tsv")

    script:
    Integer suffix = params.threshold * 100
    outputBase = [dir, "/all-tpm"].join('')
    thresholdBase = [outputBase, "t$suffix"].join("-")
    knnBase = [outputBase, '20', "k$params.knn"].join("-")
    knnThresholdBase = [outputBase, "t$suffix", "k$params.knn"].join("-")
    """
    module load MCL/$params.mclVersion

    mkdir $dir
    threshold='${params.threshold}'
    if [[ ! -z \$threshold ]]; then
        mcx alter -imx ${mci_file} -tf \
        "gq($params.threshold), add(-$params.threshold)" \
        -o ${thresholdBase}.mci
    fi
    knn='${params.knn}'
    if [[ ! -z \$knn ]]; then
        mcx alter -imx ${mci_file} -tf "add(-0.2), #knn($params.knn)" \
        -o ${knnBase}.mci
    fi

    module load datamash
    if [[ ! -z \$knn && ! -z \$threshold ]]; then
        mcx alter -imx ${mci_file} \
        -tf "gq($params.threshold), add(-$params.threshold), #knn($params.knn)" \
        -o ${knnThresholdBase}.mci
        mcx query -imx ${knnThresholdBase}.mci > ${knnThresholdBase}.stats.tsv
        # Summarise stats
        paste <( datamash -s groupby 2 count 2 < ${knnThresholdBase}.stats.tsv | awk '{ if(\$1 == 1){ print \$2 } }' ) \
            <( datamash --header-in --round 1 mean degree median degree iqr degree < ${knnThresholdBase}.stats.tsv ) | \
        awk -F"\t" 'BEGIN{ OFS = "\t" } {if(\$1 == ""){ print "0", \$2, \$3, \$4 } else{ print \$0 }}' | \
        cat <( echo -e "S\tDmean\tDmedian\tDiqr" ) - > ${knnThresholdBase}.cor-knn-stats.tsv
    fi
    """
}

// Cluster nextwork
process CLUSTER {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/all-tpm*"
    
    input:
    tuple val(dir), path(mci_file)
    each inflation

    output:
    path "*/*.mci.I[0-9]*"

    script:
    Integer inflationSuffix = inflation * 10
    """
    mkdir $dir

    module load MCL/$params.mclVersion
    mcl $mci_file -I $inflation -o $dir/${mci_file}.I${inflationSuffix}
    
    mci_base=\$( basename $mci_file .mci )
    clm info $mci_file $dir/${mci_file}.I${inflationSuffix} >> $dir/\${mci_base}.info.txt
    """
}

workflow {

    subset_output_ch = SUBSET(params.expts, params.all_counts)
        .flatten()
        .view()

    orig_ch = CREATE_BASE_NETWORK(subset_output_ch)
        .map { [it[0], it[1]] }
        .view()

    test_params_ch = TEST_PARAMETERS(orig_ch)
        .map { [it[0], it[1]] }
        .view()

    if (params.clustering) {
        threshold_ch = THRESHOLD(test_params_ch)
        .flatMap { dir = it[0]
        mci_with_id = []
        for (mci_file in it[1]) {
            mci_with_id.add([dir, mci_file])
        }
        return(mci_with_id) }
        .view()

        cluster_ch = CLUSTER(threshold_ch, params.inflationParams)
            .view()
    }
}
