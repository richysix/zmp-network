# zmp-network
Correlation network of zebrafish development

## Setup
```
scratch 
cd detct/grcz11/
export gitdir=$HOME/checkouts/zmp-network
export basedir=$(pwd)
screen -S zmp-network
```

Copy data
```
mkdir -p everything/filter-strict/
rclone copy --progress --max-depth 1 --filter "+ all.csv.gz" --filter "- *" \
sharepoint-qmul-buschlab:detct/grcz11/everything/filter-strict/ everything/filter-strict/

rclone copy --verbose --filter "+ samples.txt" --filter "- *" \
sharepoint-qmul-buschlab:detct/grcz11/everything/ everything/

# get expt and sample info
mkdir local_effect
rclone copy --verbose --filter "+ expt2samples.txt" --filter "- *" \
sharepoint-qmul-buschlab:detct/grcz11/local_effect/ local_effect/
```

Create file of expt names, sample names and conditions
```
join -t$'\t' -1 2 <( sort -t$'\t' -k2,2 local_effect/expt2samples.txt ) \
<( sort -t$'\t' -k1,1 everything/samples.txt ) | \
awk -F"\t" 'BEGIN{ OFS = "\t" } { print $2, $1, $4 }' | sort -u | \
cat <( echo -e "expt\tsample\tcondition" ) - > expt-sample-condition.tsv
# subset for testing
grep -E "^expt|(capzb1_hi1858bTg|dnmt8|zmp_ph1)\b" \
expt-sample-condition.tsv > expt-sample-condition-test.tsv
```

Create test counts file with subset of regions
```
gzip -cd everything/filter-strict/all.csv.gz | -n1001 | \
gzip -c > everything/filter-strict/all-test.csv.gz
```

## Subset by experiment

Create script to read in all counts, subset to each expt, 
aggregate counts to genes and write out subset file
```
# count lines of all file
lines=$( gzip -cd everything/filter-strict/all.csv.gz | wc -l )
python scripts/subset-by-expt.py --infer_schema_length=$lines \
expt-sample-condition-test.tsv everything/filter-strict/all-test.csv.gz
```

Create Nextflow pipeline to create and test networks for each experiment
```
mkdir nf
cd nf
# symlinks
ln -s $HOME/checkouts/zmp-network/qsub qsub
ln -s $HOME/checkouts/zmp-network/scripts scripts
ln -s $HOME/checkouts/zmp-network/nextflow.config nextflow.config
```

Run with testing set to true
```
module load nextflow
PARAMS="--testing true"
SCRIPT="scripts/main.nf"
nextflow run $PARAMS $SCRIPT
```

Run full pipeline
```
nextflow run $SCRIPT
```

## Calculate TPM
Needs transcript info. Download e98 annotation
```
cd $basedir
mkdir reference
cd reference/
wget https://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.chr.gtf.gz
```

R script to create FPKM/TPM
Expects count and sample files and transcript info
Uses the GTF to calculate transcript info and outputs a transcripts file for 
future use
e.g. for an expt dir containing samples.txt and counts-by-gene.tsv
```
Rscript scripts/counts-to-fpkm-tpm.R \
--gtf_file $basedir/reference/Danio_rerio.GRCz11.98.chr.gtf.gz \
--transcripts_file $basedir/reference/Danio_rerio.GRCz11.98.transcripts.tsv \
--fpkm --tpm --output_base $dir/all $dir/samples.txt $dir/counts-by-gene.tsv
```

## Cluster TPM

Shell script to run MCL to create and cluster a coexpression network
To run in exploratory mode to find suitable pruning parameters use -j
```
scripts/mcl-clustering-coexpr.sh -j 100/200/20 DIR TRANSCRIPT_FILE
```

To run in clustering mode use either -k (k-nearest-neighbours param) or 
-t (Correlation threshold)
```
scripts/mcl-clustering-coexpr.sh -i 1.4 -i 4 -k 150 DIR TRANSCRIPT_FILE
# or
scripts/mcl-clustering-coexpr.sh -i 1.4 -i 4 -t 0.6 DIR TRANSCRIPT_FILE
```

## Add CREATE_NETWORK process to Nextflow script

The process takes an expt directory name as input and runs 
`create-network-from-tpm.sh` with parameters set in the config file.
`create-network-from-tpm.sh` runs `counts-to-fpkm-tpm.R ` to convert the 
counts to TPM and then creates and clusters the network using 
`mcl-clustering-coexpr.sh`

Output is the same expt directory name.

## Change processes

Change CREATE_NETWORK to CREATE_BASE_NETWORK. Also change output to be the 
files created rather than the directory otherwise process doesn't get cached 
properly, because directory has being modified during each run.

Next process is TEST_PARAMETERS, which creates a base network with a 
correlation threshold of 0.2. It then runs `mcx query --vary-correlation` and
`mcx query -vary-knn`. Outputs are the mci file for the 0.2 network and the 
stats output by each of the query commands.

The next two processes are conditional on the --clustering parameter being 
set to true. This is so that the a first run of the pipeline can be done to see
the effect of varying the parameters and ten a subsequent run can be done with
threshold and clustering parameters set. 

The first process is THRESHOLD, which uses either or both of --threshold and 
-knn to reduce the edges in the network. The outputs are the mci files for each
pruned network.

The second is CLUSTERING which clusters each network for different inflation 
values set in --inflationParams. The output is the mci file for the clustered 
network.

Run the pipeline for tfap2a experiments (zmp_ph192, zmp_ph238, zmp_ph250).
```
qsub qsub/run-nextflow.sh \
-p "--expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv --knn 240 --threshold 0.44" \
scripts/main.nf 
```

Rerun with -resume and --clustering=true
```
qsub qsub/run-nextflow.sh \
-p "--expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv --knn 240 --threshold 0.44 --clustering true" \
-r current scripts/main.nf 
```
