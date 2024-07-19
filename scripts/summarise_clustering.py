#! /usr/bin/env python

import argparse
import sys
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns

def main(args):
    ''' Main body of code '''

    numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    # read in clustering file
    node_ids = []
    cluster_sizes = []
    header = True
    with open(args.input_file, mode = "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if header:
                if line[0:5] == "begin":
                    header = False
                continue
            else:
                if line[-1] == "$":
                    # cluster end
                    ids = line.split()
                    end_char = ids.pop()
                    if line[0] in numbers:
                        # cluster start
                        cluster_id = ids.pop(0)
                    node_ids.extend(ids)
                    cluster_sizes.append(len(node_ids))
                    node_ids = []
                elif line[0] in numbers:
                    # cluster start
                    ids = line.split()
                    cluster_id = ids.pop(0)
                    node_ids.extend(ids)
                else:
                    ids = line.split()
                    node_ids.extend(ids)
    # summarise clusters
    df = pl.DataFrame({"cluster_sizes": cluster_sizes})
    cl_size_fn = args.input_file + ".cl-sizes.csv"
    pl.DataFrame({
        "min": df.select(pl.col("cluster_sizes").min()),
        "min_gt1": df.filter(pl.col("cluster_sizes") > 1).select(pl.col("cluster_sizes").min()),
        "q1": df.select(pl.col("cluster_sizes").quantile(0.25)),
        "median": df.select(pl.col("cluster_sizes").median()),
        "q3": df.select(pl.col("cluster_sizes").quantile(0.75)),
        "max": df.select(pl.col("cluster_sizes").max()),
        "mean": df.select(pl.col("cluster_sizes").mean()),
        "std": df.select(pl.col("cluster_sizes").std()),
    }).write_csv(cl_size_fn, separator = "\t")

    # plot histogram of cluster sizes
    hist = sns.displot(df, x="cluster_sizes", binwidth = 10)
    cl_hist_fn = args.input_file + ".clusters.hist.png"
    hist.savefig(cl_hist_fn)

if __name__ == '__main__':
    desc = ''' Script to summarise the cluster sizes from an MCL clustering'''
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('input_file', nargs='?', metavar='INFILE',
        type=str, default='all-tpm-20-k100.mci.I14', help='Input file name')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
