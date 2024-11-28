#!/usr/bin/env bash
# subset-by-expt.sh - Script to run subset-by-expt.py
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=20G

source bash_functions.sh

USAGE="subset-by-expt.sh [options] expt_file counts_file"

OPTIONS="Options:
    -o    Output Directory [default: CWD]
    -p    Python version [default: 3.12.4]
    -s    script directory [default: $HOME/checkouts/zmp-network/scripts]
    -d    print debugging info
    -v    verbose output
    -q    turn verbose output off
    -h    print help info"

# default values
debug=0
verbose=1
OUTPUT_DIR=""
PYTHON_VERSION="3.12.4"
SCRIPT_DIR="$HOME/checkouts/zmp-network/bin"
while getopts ":o:p:s:dhqv" opt; do
  case $opt in
    o)  OUTPUT_DIR="--output_dir_prefix $OPTARG" ;;
    p)  PYTHON_VERSION=$OPTARG  ;;
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

PYTHON_VERSION=$( echo $PYTHON_VERSION | sed -e 's|^Python/||' )
module load Python/$PYTHON_VERSION

if [[ $debug -gt 0 ]]; then
    echo "Python version = $PYTHON_VERSION
Output Dir = $OUTPUT_DIR
Script Dir = $SCRIPT_DIR
Expt file: $1
Counts file: $2"
fi

if file --mime-type --brief $2 | grep -q 'gzip$'; then
  lines=$( gzip -cd $2 | wc -l )
else
  lines=$( wc -l $2 | awk '{print $1}' )
fi
subset-by-expt.py --infer_schema_length=$lines $OUTPUT_DIR $1 $2
SUCCESS=$?

verbose=1
error_checking $SUCCESS "job subset SUCCEEDED." "job nextflow FAILED: $SUCCESS"
exit $SUCCESS

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2024. Queen Mary University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
