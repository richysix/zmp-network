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

Look through stats output
There are TPM values for 26772 genes. Want threshold values with < 10% 
singletons (~2677) and median node degree of ~100.

### zmp_ph192

*tfap2a/tfap2c* incross at 4-9S.

#### Vary correlation threshold

```
column -t results/expt-zmp_ph192/all-tpm.cor-stats.tsv | \
awk '{ if(NR == 1){ print $0 } if($15 >= 0.36 && $15 <= 0.44){ print $0 }}' 
```
L      D      R     S      E           cce    EWmean    EWmed     EWiqr      NDmean   NDmed  NDiqr  CCF  eff  Cutoff
25407  1365   0     1365   0.154163    24111  0.446038  0.417025  0.0901631  486.156  235    529    NA   NA   0.36
25401  1371   0     1371   0.0948671   24100  0.48829   0.456881  0.0907996  299.164  126    317    NA   NA   0.4
25127  1645   0     1641   0.0582407   23583  0.532257  0.497171  0.0930104  183.663  68     185    NA   NA   0.44

**Good value would be between 0.4 and 0.44**

#### k-nearest neighbours

```
column -t results/expt-zmp_ph192/all-tpm.knn-stats.tsv | \
awk '{ if(NR == 1){ print $0 } if($15 >= 200 && $15 <= 280){ print $0 }}' 
```
L      D     R    S     E          cce    EWmean    EWmed     EWiqr     NDmean   NDmed  NDiqr  CCF  eff  kNN
25388  1384  0    1384  0.0439979  24075  0.483225  0.454095  0.170955  138.748  125    119    NA   NA   280
25383  1389  0    1389  0.0402976  24066  0.487921  0.459038  0.172018  127.079  113    112    NA   NA   260
25378  1394  0    1392  0.0366309  24056  0.493085  0.464496  0.173266  115.516  101    105    NA   NA   240
25368  1404  0    1400  0.0329596  24037  0.498951  0.470645  0.174534  103.939  90     97     NA   NA   220
25357  1415  0    1413  0.02933    24016  0.505518  0.477465  0.175884  92.4926  78     90     NA   NA   200

**-knn=240 looks good**

### zmp_ph238

```
column -t results/expt-zmp_ph238/all-tpm.cor-stats.tsv | \
awk '{ if(NR == 1){ print $0 } if($15 >= 0.36 && $15 <= 0.44){ print $0 }}' 
```
L      D      R     S      E            cce    EWmean    EWmed     EWiqr      NDmean    NDmed  NDiqr  CCF  eff  Cutoff
25866  906    0     906    0.149381     24990  0.463226  0.423574  0.111324   490.371   181    483    NA   NA   0.36
25853  919    0     919    0.095753     24965  0.510771  0.469036  0.124613   314.327   89     262    NA   NA   0.4
25414  1358   0     1346   0.0632589    24124  0.558243  0.516746  0.139254   207.659   43     141    NA   NA   0.44

**t=0.4 looks good**

#### k-nearest neighbours

```
column -t results/expt-zmp_ph238/all-tpm.knn-stats.tsv | \
awk '{ if(NR == 1){ print $0 } if($15 >= 200 && $15 <= 280){ print $0 }}' 
```
L      D     R    S     E          cce    EWmean    EWmed     EWiqr     NDmean   NDmed  NDiqr  CCF  eff  kNN
25687  1085  0    1060  0.040902   24646  0.446262  0.412873  0.142716  134.268  126    99     NA   NA   280
25642  1130  0    1100  0.0375998  24559  0.449998  0.416921  0.143246  123.428  115    94     NA   NA   260
25575  1197  0    1157  0.0342866  24431  0.454219  0.421404  0.143892  112.552  103    87     NA   NA   240
25502  1266  4    1219  0.0309899  24292  0.458936  0.42641   0.144805  101.73   92     80     NA   NA   220
25413  1355  4    1291  0.0277109  24123  0.464266  0.431925  0.145805  90.966   82     73     NA   NA   200

**-knn=240 looks good**

### zmp_ph250

*tfap2a/tfap2c* incross at Prim-5 (24 hpf).

#### Vary correlation threshold

```column -t results/expt-zmp_ph250/all-tpm.cor-stats.tsv | 
awk '{ if(NR == 1){ print $0 } if($15 >= 0.36 && $15 <= 0.44){ print $0 }}' 
```
L      D      R     S      E            cce    EWmean    EWmed     EWiqr        NDmean    NDmed  NDiqr  CCF  eff  Cutoff
25988  784    0     784    0.100194     25226  0.431786  0.408874  0.0783687    282.969   142    260    NA   NA   0.36
25972  800    0     800    0.0567753    25195  0.47287   0.4495    0.079552     160.346   68     148    NA   NA   0.4
25278  1494   0     1472   0.0323994    23867  0.514078  0.490278  0.0811265    91.5031   30     83     NA   NA   0.44

**-t=0.4 looks good**

#### k-nearest neighbours

```
column -t results/expt-zmp_ph250/all-tpm.knn-stats.tsv | \
awk '{ if(NR == 1){ print $0 } if($15 >= 200 && $15 <= 280){ print $0 }}' 
```
L      D    R  S    E          cce    EWmean    EWmed     EWiqr     NDmean   NDmed  NDiqr  CCF  eff  kNN
25988  784  0  784  0.0536895  25226  0.42615   0.4031    0.120983  151.631  142    83     NA   NA   280
25988  784  0  784  0.0492074  25226  0.430648  0.407289  0.122486  138.973  129    79     NA   NA   260
25988  784  0  784  0.0447728  25226  0.435545  0.411811  0.124166  126.448  116    74     NA   NA   240
25988  784  0  784  0.040361   25226  0.440924  0.416775  0.126067  113.988  104    69     NA   NA   220
25988  784  0  784  0.0359927  25226  0.446868  0.422318  0.128309  101.651  91     63     NA   NA   200

**-knn=220-240**

Rerun with -resume and --clustering=true
```
qsub qsub/run-nextflow.sh \
-p "--expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv --knn 240 --threshold 0.44 --clustering true" \
-r current scripts/main.nf 
```

```
cd ~/work/apocrita/data/scratch/bty114/detct/grcz11/nf
rsync -avzkL --filter "+ all-tpm[-.]*.[tc]sv" --filter "+ *png" --filter "+ */" \
--filter "- *" apocrita:/data/scratch/bty114/detct/grcz11/nf/results ./
```

### Testing MCL 2 nodes and edges script

Script needs mci file, gene file, cluster file and names for the output nodes
and edges files

Make annotation file with nodes ids to gene ids and names
```
# run script to download annotation
cd /data/scratch/bty114/genomes/GRCz11/e92
qsub ~/checkouts/uge-job-scripts/get_ensembl_gene_annotation.sh -e 95 \
-f ~/checkouts/bio-misc/get_ensembl_gene_annotation.pl -s "Danio rerio" 

# use annotation to create mapping from node ids to gene ids and names
join -t$'\t' -1 2 <( sort -t$'\t' -k2,2 all-tpm.tab ) \
<( awk -F"\t" 'BEGIN{OFS = "\t" } { print $1, $7, $8, NR }' Dr-e92-annotation.txt | \
sort -t$'\t' -k1,1 ) | \
awk -F"\t" 'BEGIN{OFS = "\t"} { print $2, $5, $1, $3, $4}' > node_id-gene_id-gene_name.txt
```

Download mcl2nodes script
```
cd $gitdir/scripts/
wget https://raw.githubusercontent.com/richysix/bioinf-gen/6d5ccc694b3d300fede4e5af47713dff41c111f3/mcl2nodes-edges.py
chmod a+x mcl2nodes-edges.py summarise_clustering.py
```

```
cd /data/scratch/bty114/detct/grcz11/nf-test/
mkdir tmp
cd tmp

cp ../work/4e/7015a58d55fd0fcfcac9ba857e2a2b/all-tpm-t60.mci ./
cp ../work/4e/7015a58d55fd0fcfcac9ba857e2a2b/expt-zmp_ph23/all-tpm-t60.mci.I40 ./

# create genes info file
join -t$'\t' -1 2 <( sort -t$'\t' -k2,2 all-tpm.tab ) \
<( awk -F"\t" 'BEGIN{OFS = "\t" } { print $1, $7, $8, NR }' Dr-e92-annotation.txt | \
sort -t$'\t' -k1,1 ) | awk -F"\t" 'BEGIN{OFS = "\t"} { print $2, $5, $1, $3, $4}' > node_id-gene_id-gene_name.txt

# run script
python ~/checkouts/zmp-network/scripts/mcl2nodes-edges.py \
all-tpm-t60.mci all-tpm-t60.mci.I40 node_id-gene_id-gene_name.txt \
--nodes_file expt-zmp_ph23/all-tpm-t60.I40.nodes.csv \
--edges_file expt-zmp_ph23/all-tpm-t60.I40.edges.csv --edge_offset 0.6
```

Move output to new directory
```
scratch 
mkdir zmp-network
cp -a detct/grcz11/nf zmp-network/
cd zmp-network
export gitdir=$HOME/checkouts/zmp-network
export basedir=$(pwd)
screen -S zmp-network
cd nf
```

Do fresh run of the pipeline
```
# Runtime of 240 can't be used at the moment, due to Apocrita maintenance
qsub -l h_rt=24:0:0 qsub/run-nextflow.sh \
-p "--expts $SCRATCH/zmp-network/expt-sample-condition-tfap2.tsv --knn 240 --threshold 0.44 --clustering true" \
scripts/main.nf 

# MCLTOGRAPH process did not publish the correct files
# Rerun with resume
qsub -l h_rt=24:0:0 qsub/run-nextflow.sh \
-p "--expts $SCRATCH/zmp-network/expt-sample-condition-tfap2.tsv --knn 240 --threshold 0.44 --clustering true" \
-r current scripts/main.nf 
```

Check attributes of DevStages network. Download from eLife paper
```
md5sum zfs-dev-stages-biolayout/elife-30860-supp4-v1.expression 
e1d45ce6a8f1a0539e9bc6b60bf14be8  zfs-dev-stages-biolayout/elife-30860-supp4-v1.expression
# matches
BuschLab/Documents/cam-busch/Archive/Zebrafish_Mutation_Resource/Phenotyping/Molecular%20phenotyping/RNASeq_Developmental_stages/biolayout/GRCz10/change_with_stage.expression
md5sum ~/Downloads/change_with_stage.expression 
e1d45ce6a8f1a0539e9bc6b60bf14be8  /Users/rjw26/Downloads/change_with_stage.expression
```

# Dev Stages Biolayout network

Compare to biolayout network from White et al.
Download from eLife paper
```
cd ~/work/apocrita/data/scratch/bty114/zmp-network
mkdir zfs-dev-stages-biolayout/
cd zfs-dev-stages-biolayout/

wget -O elife-30860-supp4-v1.expression \
https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzA4NjAvZWxpZmUtMzA4NjAtc3VwcDQtdjEuZXhwcmVzc2lvbg--/elife-30860-supp4-v1.expression?_hash=c5pjyRlx8sr7H9nqDJUF4H6KbC8UscObVlzYc6UjCF4%3D
```

Create script to convert biolayout file to adjancency matrix
```
# Download layout file from Sharepoint
# run script
cd ~/checkouts/zmp-network
source .venv/bin/activate
python ~/checkouts/zmp-network/scripts/biolayout-to-adj-matrix.py \
--plot_filebase ~/work/apocrita/data/scratch/bty114/zmp-network/zfs-dev-stages-biolayout/publication_change_with_stage_r-0.94 \
/Users/rjw26/work/apocrita/data/scratch/bty114/zmp-network/zfs-dev-stages-biolayout/publication_change_with_stage_r-0.94.layout \
/Users/rjw26/work/apocrita/data/scratch/bty114/zmp-network/zfs-dev-stages-biolayout/publication_change_with_stage_r-0.94.summary.txt \
/Users/rjw26/work/apocrita/data/scratch/bty114/zmp-network/zfs-dev-stages-biolayout/publication_change_with_stage_r-0.94.edge_counts.tsv 

cd ~/work/apocrita/data/scratch/bty114/zmp-network
cat zfs-dev-stages-biolayout/publication_change_with_stage_r-0.94.summary.txt
Total Edges = 1473314
Mean Edge Degree = 239.641
Quartiles
Min = 1.0
Q1 = 9.0
Median = 53.0
Q3 = 309.0
Max = 1782.0
Number of Nodes = 12296
```

## Guilt By Association

Install new package EGAD to run GBA analysis for ability of networks to predict
GO annotations (Could try ZFA as well).

Created a script `scripts/run-GBA-network.R` to run GBA in R
Added this to the MCLTOGRAPH process.

Rerun with -resume and --clustering=true
```
qsub qsub/run-nextflow.sh \
-p "--expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv --knn 240 --threshold 0.44 --clustering true -with-dag -with-report zmp-network-nf-main.html -with-timeline zmp-network-nf-timeline.html" \
-r current scripts/main.nf 
```

## Pipeline DAG

Create a mermaid diagram of the pipeline that can be included in the README.
```
qlogin
module load nextflow
nextflow run --expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv \
--knn 240 --threshold 0.44 --clustering true -preview -with-dag dag-20240923.mmd \
-resume scripts/main.nf
```

`scripts/README`

Created a Quarto document to create the README detailing the processes in the pipeline. 
Create markdown file from the mermaid diagram for including the README.
```
cat <( echo "\`\`\`mermaid" ) \
dag-20240923.mmd \
<( echo "\`\`\`" ) > $gitdir/scripts/_dag-20240923.md
```