#!/usr/bin/env nextflow

log.info """\
  NETWORK CONSTRUCTION PIPELINE
  -----------------------------

  Testing: ${params.testing}
  Expt to sample file: ${params.expts}
  All counts file: ${params.all_counts}
"""

process SUBSET {
    clusterOptions "-l h_vmem=$params.big_mem"

    input: 
    val expt_file
    val all_counts_file

    output:
    path "$params.CountDir/*"

    script:
    """ $params.QsubDir/subset-by-expt.sh -s ${params.ScriptDir} -o ${params.CountDir} ${expt_file} ${all_counts_file} """
}

process CREATE_NETWORK {
    clusterOptions "-l h_vmem=$params.big_mem"
    time { 
        if (task.attempt > 1) {
            return '240h'
        } else {
            return '1h'
        }
    }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    path expt_dir

    output:
    path expt_dir.name

    script:
    """ $params.QsubDir/create-network-from-tpm.sh \
    $params.inflationParams $params.threshold \
    $params.knnTestParams $params.knn \
    $params.outputBase $params.corMeasure $params.labels $params.skipRows \
    $params.skipCols $params.mclVersion $params.RVersion \
    $expt_dir $params.RefDir/$params.transcriptFile
    """
}




workflow {

    expt_file_ch = SUBSET(params.expts, params.all_counts)
        .flatten()
    expt_file_ch.view()

    network_ch = CREATE_NETWORK(expt_file_ch)
    network_ch.view()

}
