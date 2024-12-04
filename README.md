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
grep -E "^expt|zmp_ph(192|238|250)\b" \
expt-sample-condition.tsv > expt-sample-condition-tfap2.tsv
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
ln -s $HOME/checkouts/zmp-network/bin bin
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
-p "--expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv --knn 240 --threshold 0.44 --clustering true" \
-r current scripts/main.nf 
```

## Pipeline DAG and report

Create a mermaid diagram of the pipeline that can be included in the README.
```
qlogin
module load nextflow
cd $SCRATCH/zmp-network/nf
mkdir reports
today=$(date +%Y%m%d)
nextflow run --expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv \
--knn 240 --threshold 0.44 --clustering true -preview -with-dag reports/dag-$today.mmd \
-resume scripts/main.nf

# with report
nextflow run --expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv \
--knn 240 --threshold 0.44 --clustering true -preview -with-dag reports/dag-20240923.mmd \
-with-report reports/zmp-network-nf-main.html -with-timeline reports/zmp-network-nf-timeline.html
-resume scripts/main.nf
# if rerunning, either delete zmp-network-nf-main.html and zmp-network-nf-timeline.html
# or enable report.overwrite option in config
```

### Rerun with dag and report
```
today=$(date +%Y%m%d)
qsub qsub/run-nextflow.sh \
-p "--expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv --knn 240 --threshold 0.44 --clustering true -with-dag reports/dag-$today.mmd -with-report reports/zmp-network-nf-main-$today.html -with-timeline reports/zmp-network-nf-timeline-$today.html" \
-r current scripts/main.nf
```

`scripts/README`

Created a Quarto document to create the README detailing the processes in the pipeline. 
Create markdown file from the mermaid diagram for including the README.
```
today=$(date +%Y%m%d)
cat <( echo "\`\`\`mermaid" ) \
reports/dag-$today.mmd \
<( echo "\`\`\`" ) > $gitdir/scripts/_dag-current.md
```

### Test different thresholds effect on GO enrichment

```bash
qlogin
cd /data/scratch/bty114/zmp-network/nf/work/72/887a1c5de895b673b79f970ec0f2e6
module load MCL/14-137
mci_file=all-tpm-20.mci
dir=expt-zmp_ph192
for threshold in 0.5 0.6 0.7 0.8 0.9
do
  suffix=$( perl -le "{ print $threshold * 100 }" )
  thresholdBase="$dir/all-tpm-t$suffix"
  mcx alter -imx ${mci_file} -tf \
    "gq(${threshold}), add(-${threshold})" \
    -o ${thresholdBase}.mci
  mcx query -imx ${thresholdBase}.mci > ${thresholdBase}.stats.tsv
done

# cluster each one with lowest inflation value
inflation=1.4
inflationSuffix=14
for threshold in 0.5 0.6 0.7 0.8 0.9
do
  suffix=$( perl -le "{ print $threshold * 100 }" )
  thresholdBase="$dir/all-tpm-t$suffix"
  mci_file=${thresholdBase}.mci
  mcl $mci_file -I $inflation -o ${mci_file}.I${inflationSuffix}
done

module load Python/3.12.4
tab_file=/data/scratch/bty114/zmp-network/nf/results/expt-zmp_ph192/all-tpm.tab
annotation_file=/data/scratch/bty114/detct/grcz11/reference/Dr-e92-annotation.txt
for threshold in 0.5 0.6 0.7 0.8 0.9
do
  suffix=$( perl -le "{ print $threshold * 100 }" )
  thresholdBase="$dir/all-tpm-t$suffix"
  mci_file=${thresholdBase}.mci
  cluster_file=${mci_file}.I${inflationSuffix}
  cluster_base=$( basename $cluster_file )

  python ${gitdir}/scripts/convert_mcl.py \
--min_cluster_size 4 --graph_id ${cluster_base} \
--nodes_file ${dir}/${cluster_base}.nodes.tsv \
--edges_file ${dir}/${cluster_base}.edges.tsv \
--edge_offset ${threshold} \
$mci_file $cluster_file $tab_file $annotation_file
done

# Run GBA
qsub $gitdir/qsub/run-GBA.sh

module load datamash/1.5
for file in expt-zmp_ph192/all-tpm-t[0-9][0-9].mci.I14.go.auc.tsv
do echo -n "$file "
datamash --header-in mean 2 mean 4 < $file
zfa_file=$( echo $file | sed -e 's|go|zfa|' )
echo -n "$zfa_file "
datamash --header-in mean 2 mean 4 < $zfa_file
done | sort -k1.36,1.38
expt-zmp_ph192/all-tpm-t50.mci.I14.go.auc.tsv 0.5342999579048   0.53061049803388
expt-zmp_ph192/all-tpm-t60.mci.I14.go.auc.tsv 0.53394106524957  0.5232557758945
expt-zmp_ph192/all-tpm-t70.mci.I14.go.auc.tsv 0.53402257881198  0.51293172740254
expt-zmp_ph192/all-tpm-t80.mci.I14.go.auc.tsv 0.5363621100721   0.48506519884509
expt-zmp_ph192/all-tpm-t90.mci.I14.go.auc.tsv 0.53970784934143  0.522001939447
expt-zmp_ph192/all-tpm-t50.mci.I14.zfa.auc.tsv 0.55290938001188 0.55079021575828
expt-zmp_ph192/all-tpm-t60.mci.I14.zfa.auc.tsv 0.55108514396258 0.54324413045117
expt-zmp_ph192/all-tpm-t70.mci.I14.zfa.auc.tsv 0.56086544018161 0.53010743859729
expt-zmp_ph192/all-tpm-t80.mci.I14.zfa.auc.tsv 0.57414829893933 0.50575074120105
expt-zmp_ph192/all-tpm-t90.mci.I14.zfa.auc.tsv 0.59724296439062 0.56848400263962

cd /data/scratch/bty114/zmp-network/nf/results/expt-zmp_ph192/
for file in *.go.auc.tsv
do echo -n "$file "
datamash --header-in mean 2 mean 4 < $file
done
all-tpm-t20-k240.mci.I14.go.auc.tsv 0.52385582225004    0.4720705753021
all-tpm-t20-k240.mci.I20.go.auc.tsv 0.53278224305905    0.48778126810195
all-tpm-t20-k240.mci.I40.go.auc.tsv 0.52698166939694    0.49042056604816
all-tpm-t44-k240.mci.I14.go.auc.tsv 0.52370670419183    0.50791828242453
all-tpm-t44-k240.mci.I20.go.auc.tsv 0.52523939844065    0.507398578215
all-tpm-t44-k240.mci.I40.go.auc.tsv 0.52637475067694    0.50552012236055
all-tpm-t44.mci.I14.go.auc.tsv 0.54335982867671 0.5400352720999
all-tpm-t44.mci.I20.go.auc.tsv 0.54433722174289 0.54164264710159
all-tpm-t44.mci.I40.go.auc.tsv 0.55457281114976 0.55079479838103
```

Also try changing knn threshold
```bash
qlogin
cd /data/scratch/bty114/zmp-network/nf/work/72/887a1c5de895b673b79f970ec0f2e6
module load MCL/14-137
mci_file=all-tpm-20.mci
dir=expt-zmp_ph192
for knn in 120 140 160 180 200
do
  knnBase="$dir/all-tpm-t20-k${knn}"
  mcx alter -imx ${mci_file} -tf "add(-0.2), #knn($knn)" \
-o ${knnBase}.mci
done

# cluster each one with lowest inflation value
inflation=1.4
inflationSuffix=14
for knn in 120 140 160 180 200
do
  knnBase="$dir/all-tpm-t20-k${knn}"
  mcl $knnBase.mci -I $inflation -o $knnBase.mci.I${inflationSuffix}
done

module load Python/3.12.4
tab_file=/data/scratch/bty114/zmp-network/nf/results/expt-zmp_ph192/all-tpm.tab
annotation_file=/data/scratch/bty114/detct/grcz11/reference/Dr-e92-annotation.txt 
for knn in 120 140 160 180 200
do
  knnBase="$dir/all-tpm-t20-k${knn}"
  cluster_file=$knnBase.mci.I${inflationSuffix}
  cluster_base=$( basename $cluster_file )
  python ${gitdir}/scripts/convert_mcl.py \
--min_cluster_size 4 --graph_id ${cluster_base} \
--nodes_file ${dir}/${cluster_base}.nodes.tsv \
--edges_file ${dir}/${cluster_base}.edges.tsv \
--edge_offset 0.2 \
$knnBase.mci $cluster_file $tab_file $annotation_file
done

# Run GBA
cd expt-zmp_ph192/
qsub $gitdir/qsub/run-GBA.sh all-tpm-t20-k*.mci.I14.nodes.tsv

module load datamash/1.5
for file in  expt-zmp_ph192/all-tpm-t*k*.mci.I14.go.auc.tsv
do echo -n "$file "
datamash --header-in mean 2 mean 4 < $file
done
expt-zmp_ph192/all-tpm-t20-k120.mci.I14.go.auc.tsv 0.51677631566207     0.47159096823972
expt-zmp_ph192/all-tpm-t20-k140.mci.I14.go.auc.tsv 0.51821553010188     0.47042455024057
expt-zmp_ph192/all-tpm-t20-k160.mci.I14.go.auc.tsv 0.51950613761639     0.47010182851856
expt-zmp_ph192/all-tpm-t20-k180.mci.I14.go.auc.tsv 0.52064411389366     0.46972165940266
expt-zmp_ph192/all-tpm-t20-k200.mci.I14.go.auc.tsv 0.52106607916691     0.46993211865578
expt-zmp_ph192/all-tpm-t20-k120.mci.I14.zfa.auc.tsv 0.5205292472539     0.47113507046674
expt-zmp_ph192/all-tpm-t20-k140.mci.I14.zfa.auc.tsv 0.52243140512828    0.4704706281534
expt-zmp_ph192/all-tpm-t20-k160.mci.I14.zfa.auc.tsv 0.52492263355333    0.47036431445324
expt-zmp_ph192/all-tpm-t20-k180.mci.I14.zfa.auc.tsv 0.5262556931653     0.47022882271109
expt-zmp_ph192/all-tpm-t20-k200.mci.I14.zfa.auc.tsv 0.5274561611325     0.47040260013697
```

Run Gprofiler on clusters
```bash
$gitdir/scripts/gprofiler-on-network-clusters.R
```

Change the pipeline to make THRESHOLD 2 separate processes, FILTER_COR and FILTER_KNN
```
# thresholds as specified in nextflow.config
# params.threshold = [0.6, 0.7, 0.8, 0.9]
# params.knn = [240, 200, 160, 120, 80]
params="-with-dag -with-report -with-timeline"
options="--expts /data/scratch/bty114/detct/grcz11/expt-sample-condition-tfap2.tsv \
--clustering true"
qsub -m bea -M bty114@qmul.ac.uk qsub/run-nextflow.sh -d -o "$options" -p "$params" -r current scripts/main.nf

# or with no reports
qsub qsub/run-nextflow.sh -d -o "$options" -r current scripts/main.nf
```

```bash
module load datamash/1.5
for dir in $( find results/ -maxdepth 1 -mindepth 1 -type d )
do
  for domain in go zfa
  do
    for file in  $dir/*.$domain.auc.tsv
    do 
      echo -n $( basename $dir ) $( basename $file ) "$domain "
      datamash --header-in mean 2 mean 4 < $file | 
      awk '{ print $0, $1 - $2 }'
    done
  done
done | sed -e 's|[[:space:]]|\t|g' > auc-results.txt

for dir in $( find results/ -maxdepth 1 -mindepth 1 -type d )
do
  for domain in go zfa
  do
    for file in  $dir/*.$domain.auc.tsv
    do 
      echo -n $( basename $dir ) $( basename $file ) "$domain "
      datamash --header-in mean 2 mean 4 < $file | 
      awk '{ print $0, $1 - $2 }'
    done
  done
done | sed -e 's|[[:space:]]|\t|g' > auc-results.txt
```

Count number of clusters tested for GO enrichment
```bash
for knum in 80 120 160 200 240
do
  for inflation in 14 20 40
  do
    echo -n "$knum $inflation "
    ls -d nf/results/expt-zmp_ph192/GO/all-tpm-t20-k${knum}.mci.I${inflation}.cluster-* | \
    sed -e 's|nf/results/expt-zmp_ph192/GO/all-tpm-t20-k[0-9]*.mci.I[0-9]*.cluster-||' | sort -gr | head -n1
done
done
80 14 113
80 20 85
80 40 48
120 14 95
120 20 91
120 40 66
160 14 87
160 20 93
160 40 70
200 14 82
200 20 94
200 40 73
240 14 81
240 20 93
240 40 73

for inflation in 14 20 40
do
for knum in 80 120 160 200 240
do
    echo -n "$knum $inflation "
    ls -d nf/results/expt-zmp_ph192/GO/all-tpm-t20-k${knum}.mci.I${inflation}.cluster-* | \
    sed -e 's|nf/results/expt-zmp_ph192/GO/all-tpm-t20-k[0-9]*.mci.I[0-9]*.cluster-||' | sort -gr | head -n1
done
done
80 14 113
120 14 95
160 14 87
200 14 82
240 14 81
80 20 85
120 20 91
160 20 93
200 20 94
240 20 93
80 40 48
120 40 66
160 40 70
200 40 73
240 40 73
```

Print out definitions for columns in cor/knn stats
```
-------------------------------------------------------------------------------
 L       Percentage of nodes in the largest component
 D       Percentage of nodes in components of size at most 3 [-div option]
 R       Percentage of nodes not in L or D: 100 - L -D
 S       Percentage of nodes that are singletons
 E       Fraction of edges retained (input graph has 87883982)
 cce     Expected size of component, nodewise [ sum(sz^2) / sum^2(sz) ]
 EW      Edge weight traits (mean, median and IQR)
 ND      Node degree traits [mean, median and IQR]
 CCF     Clustering coefficient (scale 1-100)
 eff     Induced component efficiency relative to start graph (scale 1-1000)
k-NN     The knn parameter
```

Back up results directory to sharepoint (sharepoint-qmul-buschlab:Projects/zmp-network)
```bash
qsub nf/qsub/backup.sh
```

Also record the original positions of the files (symlinks)
```bash
find results/ -type l | xargs ls -lh | awk '{print $9, $11}' > results-symlinks.txt
```

## TO DO
Write a script to restore the results and work directories so that it can be
resumed without rerunning all the jobs

Test getting MCL to output a mtrix of all the correlation values
```bash
cd /data/scratch/bty114/zmp-network/nf/work/84/44dbe8bdbe5bf0d13af63aa8034ab6
mcxdump -imx expt-zmp_ph192/all-tpm-orig.mci \
-tab expt-zmp_ph192/all-tpm.tab --dump-table \
-digits 3 -sep-field "," -sep-lead "," \
-o expt-zmp_ph192/all-tpm-orig.mat.csv
```

Run script to count up cor values for histogram
```bash
# perl
cd /data/scratch/bty114/zmp-network/nf/work/84/44dbe8bdbe5bf0d13af63aa8034ab6/expt-zmp_ph192
date; perl -F"," -lane 'BEGIN{ use POSIX; %pc_counts = (); } 
{ next if $. == 1; foreach $item (@F){ 
  next if $item =~ m/^ENSDARG/; 
  $lower = floor($item*10)/10; 
  $upper = ceil($item*10)/10; 
  if($lower == $upper){ $lower = $upper - 0.1 } 
  $pc_counts{"$lower:$upper"} += 1 } } 
END{ foreach $key ( sort keys %pc_counts ){ 
  print join("\t", $key, $pc_counts{$key} ) } } ' all-tpm-orig.mat.csv | \
  sed -e 's|:|\t|' | sort -g -k1,2 > all-tpm-orig.cor-hist.txt; date 
Thu Oct 24 16:52:36 BST 2024
Thu Oct 24 17:10:33 BST 2024

#python
date; python ~/checkouts/zmp-network/scripts/cor-hist.py \
all-tpm-orig.mat.csv > all-tpm-orig.cor-hist-py.tsv ; date
Thu Oct 24 17:30:53 BST 2024
Thu Oct 24 17:41:40 BST 2024

# language 100 lines    all lines
# perl     5"           17'57"
# python   2"           10'47"
```

Select some extra experiments to test that aren't tfap2 ones
Look at expts with at least 100 sig genes and less than 550
```bash
local_effect_dir=$SCRATCH/detct/grcz11/local_effect
awk -F"\t" '{if($2 > 100 && $2 < 550){ print $0 }}' $local_effect_dir/sig_genes-expts.tsv | \
sort -k1,1 | join -t$'\t' - <( cut -f1,4 $local_effect_dir/expt_sample_stage.txt | sort -t$'\t' -u ) | \
join -t$'\t' - <( sort -t$'\t' -k1,1 $local_effect_dir/expt-info.txt ) \
> $local_effect_dir/expt-sig-count-stage-info.txt
```

Picked a few expts and looked through samples. Excluded expts with clutch
effects and row effects.
Picked 4

```
grep -E 'zmp_ph(46|71|204|213)' $local_effect_dir/expt-sig-count-stage-info.txt | \
cut -f1-3,5-8 | cat <( echo -e "expt\tsig-genes\tstage\tallele\tgene\tEnsemblID" ) - | \
column -t
expt       sig-genes  stage        allele   gene   EnsemblID
zmp_ph204  127        ZFS:0000033  hu3072   dag1   ENSDARG00000016153
zmp_ph213  464        ZFS:0000037  hu2849   neb    ENSDARG00000032630
zmp_ph46   449        ZFS:0000037  sa2042   gpaa1  ENSDARG00000074571
zmp_ph71   139        ZFS:0000037  sa18880  nod2   ENSDARG00000010756
```

Make new expt file
```
grep -E "^expt|zmp_ph(192|238|250|46|71|204|213)\b" \
$SCRATCH/detct/grcz11/expt-sample-condition.tsv > $SCRATCH/detct/grcz11/expt-sample-condition-tfap2-plus.tsv
```

Count samples per expt
```
cut -f1 $SCRATCH/detct/grcz11/expt-sample-condition-tfap2-plus.tsv | \
grep -v expt | uniq -c
     90 zmp_ph192
     23 zmp_ph204
     24 zmp_ph213
     87 zmp_ph238
     89 zmp_ph250
     24 zmp_ph46
     22 zmp_ph71
```

Run new pipeline with email notfication
```
qsub -m bea -M bty114@qmul.ac.uk qsub/run-nextflow.sh -d \
-o "--expts $SCRATCH/detct/grcz11/expt-sample-condition-tfap2-plus.tsv --clustering true" \
-p "-with-dag -with-report -with-timeline" scripts/main.nf
# to rerun
qsub -m bea -M bty114@qmul.ac.uk qsub/run-nextflow.sh -d \
-o "--expts $SCRATCH/detct/grcz11/expt-sample-condition-tfap2-plus.tsv --clustering true" \
-p "-with-dag -with-report -with-timeline" -r current scripts/main.nf

# or with no reports
qsub qsub/run-nextflow.sh -d -o "$options" -r current scripts/main.nf
```

Write script to get GO annotation
```
# Download script from uge repo
cd $basedir/nf/bin
wget https://github.com/richysix/uge-job-scripts/raw/e5f5faad28b23ff55419726beef3675f9e5fdba3/get-go-annotation.sh
chmod a+x get-go-annotation.sh 
cd $basedir/nf/
mkdir reference
cd reference/
../bin/get-go-annotation.sh -e 109 -s danio_rerio
```

Update config to use profiles and module/process specific options
```
mkdir conf
```

Download base nf-core config from Ian's zfvarcall repo as start and edit
And apocrita config
```
cd conf
wget https://github.com/iansealy/zfvarcall/raw/refs/heads/master/conf/base.config
wget https://github.com/iansealy/zfvarcall/raw/refs/heads/master/conf/apocrita.config
```
