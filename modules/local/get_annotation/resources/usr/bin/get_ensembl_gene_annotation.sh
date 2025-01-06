#!/usr/bin/env bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -o get_ensembl_gene_annotation.o
#$ -e get_ensembl_gene_annotation.e

USAGE="get_ensembl_gene_annotation.sh [options]"

source bash_functions.sh

OPTIONS="Options:
    -e    Ensembl version number [109]
    -f    path to get annotation script [get_ensembl_gene_annotation.pl]
    -s    Species name [Danio rerio]
    -o    Output file name [annotation.txt]
    -d    print debugging info
    -v    verbose output [default]
    -q    quiet output
    -h    print help info"

# default values
debug=0
verbose=1
ENSEMBL_VERSION=109
SPECIES='Danio rerio'
OUTPUT_FILE='annotation.txt'
FILE=get_ensembl_gene_annotation.pl

while getopts ":e:f:s:o:dvqh" opt; do
  case $opt in
    e) ENSEMBL_VERSION=$OPTARG ;;
    f) FILE=$OPTARG ;;
    s) SPECIES=$OPTARG ;;
    o) OUTPUT_FILE=$OPTARG ;;
    d) debug=1 ;;
    v) verbose=1 ;;
    q) verbose=0 ;;
    h)  usage "" ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

if [[ $debug -gt 0 ]];then
    echo "Version = $ENSEMBL_VERSION" >&2
    echo "Script File = $FILE" >&2
    echo "Species = $SPECIES" >&2
    echo "Output file = $OUTPUT_FILE" >&2
fi

module purge
module load Ensembl/$ENSEMBL_VERSION

perl $FILE --species $SPECIES > $OUTPUT_FILE

SUCCESS=$?

error_checking $SUCCESS "job get_ensembl_gene_annotation SUCCEEDED." "job get_ensembl_gene_annotation FAILED: $SUCCESS"
exit $SUCCESS
