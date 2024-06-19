#!/usr/bin/env bash
# create-network-from-tpm.sh - Script to create TPM and
# then use MCL to create a basic correlation network
# Parameters are 
# Directory: This should contain counts-by-gene.tsv and samples.tsv
# Transcripts_file: A file containing GeneIDs and gene lengths
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=20G

source bash_functions.sh

USAGE="create-network-from-tpm.sh [options] Directory Transcripts_file"

OPTIONS="Options:
    -i    Inflation values. Can be specified multiple times [default: 1.4]
    -j    Testing effect of k-nearest-neighbours. 
          Str "START/END/STEP" [default: 50/500/50]
    -k    Value for k-nearest-neighbours. [default: None]
    -o    Output basename [default: input file name minus .tsv]
    -t    Edge weight threshold [default: 0]
          If neither -k or -t is set, -t=0.6 will be used
    -c    Coexpression measure [pearson, spearman, cosine]
    -l    Column of labels [default: 1]
    -r    Skip first n rows [default: 1]
    -s    Skip first n columns [default: 1]
    -m    MCL version [default: 14-137]
    -v    R version [default: 4.3.1]
    -d    script directory [default: $HOME/checkouts/zmp-network/scripts]
    -h    print help info"

# default values
INFL_PARAMS=()
THRESHOLD=0
MEASURE="pearson"
LABELS="1"
SKIP_ROWS="1"
SKIP_COLS="1"
KNN_PARAMS="50/500/50"
KNN=""
MCL_VERSION="14-137"
R_VERSION="4.3.1"
SCRIPT_DIR="$HOME/checkouts/zmp-network/scripts"
while getopts ":i:j:k:o:t:c:l:r:s:m:v:d:h" opt; do
  case $opt in
    i)  INFL_PARAMS+=("$OPTARG") ;;
    j)  KNN_PARAMS=$OPTARG ;;
    k)  KNN=$OPTARG ;;
    o)  OUTPUT_BASE=$OPTARG  ;;
    t)  THRESHOLD=$OPTARG  ;;
    c)  MEASURE=$OPTARG ;;
    l)  LABELS=$OPTARG ;;
    r)  SKIP_ROWS=$OPTARG ;;
    s)  SKIP_COLS=$OPTARG ;;
    m)  MCL_VERSION=$OPTARG ;;
    v)  R_VERSION=$OPTARG  ;;
    d)  SCRIPT_DIR=$OPTARG ;;
    h)  usage "" ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

# check options
if [[ ! -z $KNN ]] && [[ $THRESHOLD != 0 ]]; then
  echo "Both -k and -t are set. Please set one or the other" >&2
  exit 2
elif [[ ! -z $KNN ]]; then
  FILTER_OPT="-k $KNN"
elif [[ $THRESHOLD != 0 ]]; then
  FILTER_OPT="-t $THRESHOLD"
else
  echo "Neither -k or -t is set. -t=0.6 will be used" >&2
  FILTER_OPT="-t 0.6"
fi

R_VERSION=$( echo $R_VERSION | sed -e 's|^R/||' )
module load R/$R_VERSION

# unpack arguments
dir=$1
transcript_file=$2

# create samples file
awk -F"\t" '{if(NR > 1){ print $2 "\t" $3 }}' \
 $dir/samples.tsv > $dir/samples.txt

# run counts-to-fpkm-tpm
Rscript $SCRIPT_DIR/counts-to-fpkm-tpm.R \
--transcripts_file $transcript_file \
--output_base $dir/all --output_format tsv \
--tpm $dir/samples.txt $dir/counts-by-gene.tsv

SUCCESS=$?
error_checking $SUCCESS "counts-to-tpm SUCCEEDED." "counts-to-tpm FAILED: $SUCCESS"

# run inital MCL
$SCRIPT_DIR/mcl-clustering-coexpr.sh \
-j $KNN_PARAMS -c $MEASURE -l $LABELS -r $SKIP_ROWS -s $SKIP_COLS -m $MCL_VERSION \
$dir/all-tpm.tsv

# run MCL clustering
INFLATION=""
for INFL in "${INFL_PARAMS[@]}"; do
  INFLATION="$INFLATION -i $INFL"
done
$SCRIPT_DIR/mcl-clustering-coexpr.sh \
$INFLATION $FILTER_OPT -c $MEASURE -l $LABELS -r $SKIP_ROWS -s $SKIP_COLS -m $MCL_VERSION \
$dir/all-tpm.tsv
