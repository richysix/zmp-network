#!/usr/bin/env bash
# get-go-annotation.sh - Script to get GO annotation file
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G

source bash_functions.sh

USAGE="get-go-annotation.sh [options]"

OPTIONS="Options:
    -e    Ensembl version number [109]
    -s    Species name [danio_rerio]
    -d    Turn debugging on
    -q    quiet output
    -h    print help info"

# default values
debug=0
verbose=1
ENSEMBL_VERSION=109
SPECIES='danio_rerio'
while getopts ":e:s:dqh" opt; do
  case $opt in
    e) ENSEMBL_VERSION=$OPTARG ;;
    s) SPECIES=$OPTARG ;;
    d) debug=1 ;;
    q) verbose=0 ;;
    h)  usage "" ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

BASE_URL="https://raw.githubusercontent.com/iansealy/topgo-wrapper/refs/heads/master/data"
FILE="${SPECIES}_e${ENSEMBL_VERSION}_go.txt"
GO_URL="$BASE_URL/$FILE"

if [[ $verbose -eq 1 ]]; then
  CMD="wget $GO_URL"
else
  CMD="wget -q $GO_URL"
fi
if [[ $debug -gt 0 ]]; then
  echo $CMD >&2
fi
eval $CMD
SUCCESS=$?
PATH="$(pwd)/$FILE"
error_checking $SUCCESS "GO file: $PATH" "Downloading GO file FAILED: $SUCCESS"

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
