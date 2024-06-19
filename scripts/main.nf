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

workflow {

    expt_file_ch = SUBSET(params.expts, params.all_counts)
        .flatten()
    expt_file_ch.view()

    expts_ch = expt_file_ch
        .splitText()
        .map { it.trim() }
        .view()

}
