#!/usr/bin/env bash
# run-GBA.sh - Script to run Guil-By-Association analysis
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=240:0:0
#$ -l h_vmem=128G

source bash_functions.sh

USAGE="run-GBA.sh [options] node_files"

OPTIONS="Options:
    -g    GO annotation file [default: $SCRATCH/detct/grcz11/reference/danio_rerio_e85_go.txt]
    -z    ZFA annotation file [default: $SCRATCH/detct/grcz11/reference/Dr-e98-Gene2ZFA.txt]
    -r    R version [default: 4.3.1]
    -s    script directory [default: $HOME/checkouts/zmp-network/scripts]
    -d    print debugging info
    -v    verbose output
    -q    turn verbose output off
    -h    print help info"

# default values
debug=0
verbose=1
R_VERSION="4.3.1"
SCRIPT_DIR="$HOME/checkouts/zmp-network/scripts"
go_annotation_file=$SCRATCH/detct/grcz11/reference/danio_rerio_e85_go.txt
zfa_annotation_file=$SCRATCH/detct/grcz11/reference/Dr-e98-Gene2ZFA.txt
while getopts ":g:z:r:s:dhqv" opt; do
  case $opt in
    g)  go_annotation_file=$OPTARG  ;;
    z)  zfa_annotation_file=$OPTARG  ;;
    r)  R_VERSION=$OPTARG  ;;
    s)  SCRIPT_DIR=$OPTARG ;;
    d)  debug=1  ;;
    h)  usage "" ;;
    v)  verbose=1 ;;
    q)  verbose=0 ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

module load R/$R_VERSION

inflation=1.4
inflationSuffix=14

for file in $@
do
  cluster_base=$( basename -s .nodes.tsv $file )

  if [[ $debug -gt 0 ]]; then
    echo "Progam args: $@"
    echo "First file: $file"
    echo "Basename: $cluster_base"
  fi

  Rscript /data/home/bty114/checkouts/zmp-network/scripts/run-GBA-network.R \
--auc_file ${cluster_base}.go.auc.tsv \
--scores_file ${cluster_base}.go.gene-scores.tsv \
--plots_file ${cluster_base}.go.GBA-plots.pdf \
$file \
${cluster_base}.edges.tsv \
$go_annotation_file

  Rscript /data/home/bty114/checkouts/zmp-network/scripts/run-GBA-network.R \
--auc_file ${cluster_base}.zfa.auc.tsv \
--scores_file ${cluster_base}.zfa.gene-scores.tsv \
--plots_file ${cluster_base}.zfa.GBA-plots.pdf \
$file \
${cluster_base}.edges.tsv \
$zfa_annotation_file

done
