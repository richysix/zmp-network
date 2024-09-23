# ZMP Network

## Processes

``` mermaid
flowchart TB
    subgraph " "
    v0["expt_file"]
    v1["all_counts_file"]
    v15["infl_values_ch"]
    v21["annotation_ch"]
    v22["go_annotation_ch"]
    end
    v2([SUBSET])
    v5([CREATE_BASE_NETWORK])
    v9([TEST_PARAMETERS])
    v12([THRESHOLD])
    v16([CLUSTER])
    v23([MCLTOGRAPH])
    subgraph " "
    v25[" "]
    end
    v3(( ))
    v6(( ))
    v10(( ))
    v13(( ))
    v24(( ))
    v0 --> v2
    v1 --> v2
    v2 --> v3
    v3 --> v5
    v5 --> v6
    v6 --> v9
    v9 --> v10
    v10 --> v12
    v12 --> v13
    v15 --> v16
    v13 --> v16
    v16 --> v6
    v21 --> v23
    v22 --> v23
    v6 --> v23
    v23 --> v24
    v24 --> v25
```

### SUBSET (big_mem_retry)

Creates a sample and count file for each individual experiment with the
counts aggregated to the gene level.

Inputs:

1.  Experiment file: Tab-separated file, containing columns `expt` and
    `sample`

2.  Counts file: Comma-separated file containing counts for all the
    samples (Can be gzipped)

Outputs:

1.  Tuple of directory names for each experiment with a DirPrefix (from
    the config file) added.

This gets flattened so that the next process gets the directory names
one at a time and runs `CREATE_BASE_NETWORK` once for each dir name.

Script: Runs `scripts/subset-by-expt.py`

### CREATE_BASE_NETWORK

Creates a base correlation network containing ALL edges.

Inputs:

1.  Expt dir: path to experiment directory. Contains subset sample and
    count files for the experiment.

Outputs:

1.  Dir name: name of experiment

2.  TPM file: Tab-separated file of TPMs used to create the network

3.  TAB file: File output by MCL mapping node ids to gene ids

4.  MCI file: MCL network file detailing all the edges in the network
    and their weights.

Script: Runs `scripts/counts-to-fpkm-tpm.R` and then `mcxarray`

### TEST_PARAMETERS

Creates a basic network thresholded at 0.2 and then collects stats on
node degree and singletons when the threshold is varied. Also collects
stats on using k-nearest neighbours.

Inputs:

1.  Dir name: name of experiment

2.  TPM file: Tab-separated file of TPMs used to create the network

Outputs:

1.  Dir name: name of experiment

2.  MCI file: MCL network file detailing all the edges in the network
    and their weights.

3.  Cor stats file: File of node stats varying correlation threshold

4.  Knn stats file: File of node stats varying knn threshold

Script: Runs `mcxarray` and `mcx query`

### THRESHOLD

Using the threshold and knn parameters to create pruned networks. If
both correlation and knn thresholds are set, a third network using both
thresholds together is created.

Inputs:

1.  Dir name: name of experiment

2.  MCI file: MCL network file detailing all the edges in the network
    and their

Outputs:

1.  Dir name: name of experiment

2.  Tuple of MCI files: MCL network files pruned using correlation
    and/or knn threshold. These are flattened so that the next process
    receives the expt name, string of threshold parameters and mci file
    for those parameters.

3.  Stats file created if both correlation and knn threshold are set

Script: Runs `mcx alter` and `mcx query`

### CLUSTER

Clusters the supplied network.

    tuple val(dir), val(threshold), path(mci_file)
    each inflation

Inputs:

1.  Tuple of

    1.  Dir name: name of experiment

    2.  Threshold string

    3.  MCI file

2.  Inflation values: From `inflationParams` in the config. This is a
    list. The CLUSTER process is run once for each inflation value in
    the list.

Outputs:

Tuple of

1.  Dir name: name of experiment

2.  Threshold value

3.  Inflation value

4.  MCI file

5.  Clustered MCI file

6.  Stats file: info on each node

7.  Tuple of summarise_clustering.py outputs (cluster sizes and
    histograms)

The outputs are crossed with the output from CREATE_BASE_NETWORK and
rearranged so that the original tab file is available for the next
process (MCLTOGRAPH)

Script: Runs `mcl`, `clm info` and `scripts/summarise_clustering.py`

### MCLTOGRAPH

Converts the clustered MCI file to nodes and edges files. Then runs GBA
on the network.

Inputs:

1.  Tuple of:

    1.  Dir name: name of experiment

    2.  TAB file

    3.  Threshold value

    4.  Inflation value

    5.  MCI file

    6.  Clustered MCI file

2.  Annotation file: Tab-separated file of Gene annotation (Chr, Start
    End, ID, Name etc.)

3.  GO annotation file: Tab-separated file of GO annotation (GeneID,
    TermID, Component)

Outputs:

Tuple of

1.  Dir name: name of experiment

2.  Threshold value

3.  Inflation value

4.  Nodes file: Comma-separated file of node information with Name and
    Cluster id

5.  Edges files: Comma-separated file of edge information with source,
    target and weight

6.  AUC file: Output file from GBA with AUC values for each GO term

7.  Gene scores file: Output file from GBA with gene scores for each GO
    term

8.  Plots file: Histograms of AUC values

Script: Runs `scripts/mcl2nodes-edges.py` and
`scripts/run-GBA-network.R`
