#!/usr/bin/env bash
# mcl-clustering-coexpr.sh - Script to create a clustered network from RNA-seq
# expression levels
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=12G
#$ -o mcl-clustering.o
#$ -e mcl-clustering.e

source bash_functions.sh

USAGE="Script to create a clustered correlation network

mcl-clustering-coexpr.sh [options] input_file

input_file: tab-separated file with expression counts
"

OPTIONS="Options:
    -i    Inflation values. Can be specified multiple times [default: 1.4]
    -j    Testing effect of k-nearest-neighbours. 
          Str "START/END/STEP" [default: 50/500/50]
    -k    Value for k-nearest-neighbours. [default: None]
    -o    Output basename [default: input file name minus .tsv]
    -t    Edge weight threshold [default: 0]
    -c    Coexpression measure [pearson, spearman, cosine]
    -l    Column of labels [default: 1]
    -r    Skip first n rows [default: 1]
    -s    Skip first n columns [default: 1]
    -m    MCL version [default: 14-137]
    -d    print debugging info
    -v    verbose output
    -q    turn verbose output off
    -h    print help info"

# default values
debug=0
verbose=1
INFL_PARAMS=()
THRESHOLD=0
MEASURE="--pearson"
LABELS="-l 1"
SKIP_ROWS="-skipr 1"
SKIP_COLS="-skipc 1"
KNN_PARAMS="50/500/50"
KNN=""
MCL_VERSION="14-137"

while getopts ":i:j:k:o:t:c:m:l:r:s:dhqv" opt; do
  case $opt in
    i)  INFL_PARAMS+=("$OPTARG") ;;
    j)  KNN_PARAMS=$OPTARG ;;
    k)  KNN=$OPTARG ;;
    o)  OUTPUT_BASE=$OPTARG  ;;
    t)  THRESHOLD=$OPTARG  ;;
    c)  MEASURE="--$OPTARG" ;;
    m)  MCL_VERSION="$OPTARG" ;;
    l)  LABELS="-l $OPTARG" ;;
    r)  SKIP_ROWS="-skipr $OPTARG" ;;
    s)  SKIP_COLS="-skipc $OPTARG" ;;
    d)  debug=1  ;;
    h)  usage "" ;;
    v)  verbose=1 ;;
    q)  verbose=0 ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

# unpack arguments
if [[ ! -z $KNN ]] && [[ $THRESHOLD != 0 ]]; then
  echo "Both -k and -t are set. Please set one or the other" >&2
  exit 2
fi

# define output base it it's empty
INPUT_FILE=$1
if [[ -z $OUTPUT_BASE ]]; then
    OUTPUT_BASE=$( echo $INPUT_FILE | sed -e 's|\.tsv||' )
fi

# set default inflation params if empty
if [[ ${#INFL_PARAMS[@]} -eq 0 ]]; then
  INFL_PARAMS=(1.4)
fi
# load MCL module
MCL_VERSION=$( echo $MCL_VERSION | sed -e 's|^MCL/||' )
module load MCL/$MCL_VERSION

# CLUSTER FUNCTION
cluster()
{
  cat /dev/null > $MCI_BASE.info.txt
  for INFL in "${INFL_PARAMS[@]}"
  do
    INFLATION_SUFFIX=$( echo $INFL | sed -e 's|\.|-|' )
    mcl $MCI_BASE.mci -I $INFL -o $MCI_BASE.mci.I${INFLATION_SUFFIX}

    clm info $MCI_BASE.mci \
    $MCI_BASE.mci.I${INFLATION_SUFFIX} >> $MCI_BASE.info.txt
  done

  clm dist --chain $MCI_BASE.mci.I* > $MCI_BASE.dist.txt

  SUCCESS=$?
}

# create filtered network
if [[ "$THRESHOLD" != 0 ]]; then
  # prune edges below the threshold
  SUFFIX=$( printf '%.0f\n' $( echo "$THRESHOLD * 100" | bc -l ) )
  MCI_BASE="${OUTPUT_BASE}-${SUFFIX}"
  mcx alter -imx ${OUTPUT_BASE}-20.mci -tf "gq($THRESHOLD), add(-$THRESHOLD)" \
  -o $MCI_BASE.mci
  # run cluster function
  cluster
  SUCCESS=$?
elif [[ ! -z $KNN ]]; then
  # run k-nearest neighbours
  SUFFIX="k$KNN"
  MCI_BASE="${OUTPUT_BASE}-20-${SUFFIX}"
  mcx alter -imx ${OUTPUT_BASE}-20.mci -tf "add(-0.2), #knn($KNN)" \
  -o $MCI_BASE.mci
  # run cluster function
  cluster
  SUCCESS=$?
else
  # create base network, no adjustments
  mcxarray -data $INPUT_FILE -co 0 $SKIP_ROWS $SKIP_COLS \
  $MEASURE $LABELS -o ${OUTPUT_BASE}-orig.mci -write-tab $OUTPUT_BASE.tab

  mcxarray -data $INPUT_FILE -co 0.2 $SKIP_ROWS $SKIP_COLS -tf 'abs()' \
  $MEASURE $LABELS -o ${OUTPUT_BASE}-20.mci
  
  # vary correlation
  mcx query -imx ${OUTPUT_BASE}-20.mci --vary-correlation --output-table \
   > ${OUTPUT_BASE}.cor-stats.tsv

  # test varying k-nearest neighbours
  mcx query -imx ${OUTPUT_BASE}-20.mci -vary-knn $KNN_PARAMS --output-table \
   > ${OUTPUT_BASE}.knn-stats.tsv
  # finish
  SUCCESS=$?
fi
  
error_checking $SUCCESS "job mcl SUCCEEDED." "mcl FAILED: $SUCCESS"
exit $SUCCESS

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2023. Queen Mary University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
