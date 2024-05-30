#! /usr/bin/env python

''' Script to take a DETCT counts file, subset to required samples
and aggregate the count data from regions to genes '''

import argparse
import sys
import polars.selectors as cs
import polars as pl
import re
import gzip

def main(args):
    ''' Main body of code '''

    # read in samples file
    if args.subset_samples is not None:
        sample_names = pl.read_csv(args.subset_samples, separator="\t"
            ).select("sample").with_columns(pl.Series("suffix", [" count"])
            ).select(pl.concat_str(pl.all())
            ).to_series(0).to_list()

    # read in input
    if re.search("\\.gz$", args.input_file):
        fh = gzip.open(args.input_file)
    else:
        fh = open(args.input_file)

    df = pl.read_csv(fh, infer_schema_length = args.infer_schema_length)
    if args.subset_samples is not None:
        df = df.select(cs.by_name(df.columns[0:15]) | 
                cs.by_name(sample_names)
            )
    else:
        df = df.select(~cs.ends_with("normalised count"))

    col_map = {}
    for col_name in df.columns:
        if re.match("e[0-9]+ Ensembl Gene ID", col_name):
            col_map[col_name] = "GeneID"

    if col_map:
        df = df.rename(col_map)

    # remove rows where Gene ID is null or there are multiple IDs
    # group by gene and sum counts
    df = df.filter(pl.col("GeneID").str.contains("ENSDARG")
        ).filter(~pl.col("GeneID").str.contains(",")
        ).select(~cs.contains("3' end read count")
        ).group_by(pl.col("GeneID")
        ).agg(
            cs.ends_with("count").sum()
        )
    df.write_csv(args.output_file)

if __name__ == '__main__':
    desc = ''' Script to take a DETCT counts file, subset to required samples
and aggregate the count data from regions to genes '''
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('input_file', nargs='?', metavar='INFILE',
        type=str, default='all.csv', help='Input file name')
    parser.add_argument('output_file', nargs='?', metavar='OUTFILE',
        type=argparse.FileType('w'), default=sys.stdout, help='Output file name')
    parser.add_argument('--subset_samples', default=None,
        help='File name of samples to subset to')
    parser.add_argument('--infer_schema_length', default=1000, type=int,
        help='Alter the number of rows to use to infer the data schema')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
