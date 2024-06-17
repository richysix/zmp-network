#!/usr/bin/env nextflow

params.testing = false
params.BaseDir = "/data/scratch/bty114/detct/grcz11"
if (params.testing) {
    params.expts = "$params.BaseDir/expt-sample-test.txt"
    params.all_counts = "$params.BaseDir/everything/filter-strict/all-test.csv.gz"
    params.big_mem = '1G'
} else {
    params.expts = "$params.BaseDir/expt-sample.txt"
    params.all_counts = "$params.BaseDir/everything/filter-strict/all.csv.gz"
    params.big_mem = '16G'
}
params.ScriptDir = "/data/home/bty114/checkouts/zmp-network/scripts"
params.QsubDir = "/data/home/bty114/checkouts/zmp-network/qsub"

log.info """\
  NETWORK CONSTRUCTION PIPELINE
  -----------------------------

  Testing: ${params.testing}
  Expt to sample file: ${params.expts}
  All counts file: ${params.all_counts}
"""

process SUBSET {
    clusterOptions "-l h_vmem=$params.big_mem"
    cpus 1
    time '1h'
    penv 'smp'

    input: 
    val expt_file
    val all_counts_file

    output:
    path 'expts.txt'

    script:
    """ $params.QsubDir/subset-by-expt.sh -s ${params.ScriptDir} ${expt_file} ${all_counts_file} """
}

workflow {

    expt_file_ch = SUBSET(params.expts, params.all_counts)
    expt_file_ch.view()

    expts_ch = expt_file_ch
        .splitText()
        .map { it.trim() }
        .view()

}
