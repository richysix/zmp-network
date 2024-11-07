#! /usr/bin/env python

''' Script to take a DETCT counts file, subset to required samples
and aggregate the count data from regions to genes.
It expects a samples file with experiment names and sample names and
a count file 
This script creates directories for each experiment in the current working
directory if they don't already exist and saves a file of counts, subset
to the correct samples and aggregated to the gene level. '''

import argparse
import polars.selectors as cs
import polars as pl
import re
import gzip
import os
import sys
from itertools import repeat

def output_aggregate_counts_for_expt(expt_name, df, sample_info, expts_outfh, args):
    ''' Function to subset the data frame to just the counts for an 
    experiment aggregate the counts from the transcript to the gene 
    level and output to a file '''
    # get samples for expt
    sample_names = sample_info.filter(expt = expt_name
        ).select("sample").with_columns(pl.Series("suffix", [" count"])
        ).select(pl.concat_str(pl.all())
        ).to_series(0).to_list()

    # check if any of the samples are present in the count file
    if not any([ x in df.columns for x in sample_names ]):
        print(f"{expt_name}: None of the required samples are present in "
                f"the counts file", file = sys.stderr)
        return(None)

    # create directory
    if args.output_dir_prefix is None:
        output_dir = expt_name
    else:
        output_dir = "-".join([args.output_dir_prefix, expt_name])
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    print(f"{output_dir}", file = expts_outfh)
    
    # write samples file
    outfile = os.path.join(output_dir, "samples.tsv")
    sample_info.filter(expt = expt_name
        ).write_csv(outfile, separator = "\t")
        
    # select columns by name
    subset = df.select(cs.by_name(df.columns[0:15]) | 
            cs.by_name(sample_names)
        ).filter(pl.col("GeneID").str.contains("ENSDARG") & 
        ~pl.col("GeneID").str.contains(",")
        ).select(~cs.contains("3' end read count")
        ).group_by(pl.col("GeneID", "Gene name")
        ).agg(
            cs.ends_with("count").sum()
        )
    # output counts file
    outfile = os.path.join(output_dir, "counts-by-gene.tsv")
    subset.write_csv(outfile, separator = "\t")

def main(args):
    ''' Main body of code '''

    # read in samples file
    sample_info = pl.read_csv(args.sample_file, separator="\t")
    expts = sample_info.unique(subset = "expt"
        ).select("expt"
        ).to_series(0).to_list()

    # read in input
    if re.search("\\.gz$", args.count_file):
        fh = gzip.open(args.count_file)
    else:
        fh = open(args.count_file)

    all_data = pl.read_csv(fh, infer_schema_length = args.infer_schema_length)

    col_map = {}
    for col_name in all_data.columns:
        if re.match("e[0-9]+ Ensembl Gene ID", col_name):
            col_map[col_name] = "GeneID"

    if col_map:
        all_data = all_data.rename(col_map)

    # open output file
    expts_outfh = open(args.expts_outfile, mode = 'w')

    for expt in expts:
        if args.verbose:
            print(expt)
        output_aggregate_counts_for_expt(expt, all_data, sample_info, expts_outfh, args)

if __name__ == '__main__':
    desc = ''' Script to take a DETCT counts file, subset to required samples
    and aggregate the count data from regions to genes '''
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('sample_file', nargs='?', metavar='SAMPLES',
        type=str, default='samples.tsv', help='Samples file name')
    parser.add_argument('count_file', nargs='?', metavar='COUNTS',
        type=str, default='all.csv', help='Counts file name')
    parser.add_argument('--output_dir_prefix', default=None, type=str,
        metavar = "OUTPUT DIR PREFIX",
        help='Directory to output expts to [default: %(default)s]')
    parser.add_argument('--expts_outfile', default='expts.txt', type=str,
        metavar = "EXPTS FILE",
        help='File name to output expt names to [default: %(default)s]')
    parser.add_argument('--infer_schema_length', default=1000, type=int,
        metavar = "ROWS",
        help=''.join(['Alter the number of rows to use to infer the data ',
        'schema [default: %(default)s]']))
    parser.add_argument('--verbose', action='store_true', default=False,
        help='Turns on verbose output')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
