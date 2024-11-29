process FILTER_COR {
    label 'big_mem_retry'
    publishDir "results", pattern: "*/*all-tpm-t*"
    
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

// parameters for testing
mci_file = file("$params.repo/data/all-tpm-orig.mcx")
dir = 'test-1'
workflow {
    FILTER_COR(Channel.of([dir, mci_file]), Channel.of(0.4))
    FILTER_COR.out.view()
}