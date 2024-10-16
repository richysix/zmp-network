#! /usr/bin/env python

import argparse
import sys
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

def main(args):
    ''' Main body of code '''

    numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    # read in clustering file
    node_ids = []
    cluster_ids = []
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
                        cluster_ids.append(cluster_id)
                    node_ids.extend(ids)
                    cluster_sizes.append(len(node_ids))
                    node_ids = []
                elif line[0] in numbers:
                    # cluster start
                    ids = line.split()
                    cluster_id = ids.pop(0)
                    cluster_ids.append(cluster_id)
                    node_ids.extend(ids)
                else:
                    ids = line.split()
                    node_ids.extend(ids)
    # summarise clusters
    df = pl.DataFrame({
        "cluster_id": cluster_ids,
        "cluster_size": cluster_sizes 
    })
    cl_size_fn = args.input_file + ".cl-sizes.tsv"
    df.write_csv(cl_size_fn, separator = "\t")
    no_singletons = df.filter(pl.col("cluster_size") > 1)

    cl_size_fn = args.input_file + ".cl-size-summary.tsv"
    pl.DataFrame({
        "min": df.select(pl.col("cluster_size").min()),
        "q1": df.select(pl.col("cluster_size").quantile(0.25)),
        "median": df.select(pl.col("cluster_size").median()),
        "q3": df.select(pl.col("cluster_size").quantile(0.75)),
        "max": df.select(pl.col("cluster_size").max()),
        "mean": df.select(pl.col("cluster_size").mean().round(2)),
        "std": df.select(pl.col("cluster_size").std().round(2)),
        "min_gt1": no_singletons.select(pl.col("cluster_size").min()),
        "q1_gt1": no_singletons.select(pl.col("cluster_size").quantile(0.25)),
        "median_gt1": no_singletons.select(pl.col("cluster_size").median()),
        "q3_gt1": no_singletons.select(pl.col("cluster_size").quantile(0.75)),
        "max_gt1": no_singletons.select(pl.col("cluster_size").max()),
        "mean_gt1": no_singletons.select(pl.col("cluster_size").mean().round(2)),
        "std_gt1": no_singletons.select(pl.col("cluster_size").std().round(2)),
    }).write_csv(cl_size_fn, separator = "\t")

    # plot histogram of cluster sizes
    cl_hist_fn = args.input_file + ".clusters.hist.png"
    try:
        hist = sns.displot(df, x="cluster_size", binwidth = 10)
        hist.savefig(cl_hist_fn)
    except ValueError:
        print("There was a problem with the cluster size",
            f"histogram, {cl_hist_fn}", file = sys.stderr)
        

    # plot histogram with singletons removed
    # displot throws an error if there is only one value
    if no_singletons.height <= 1:
        print("There aren't enough clusters larger than 1 to plot a histogram",
            file = sys.stderr)
    else:
        cl_hist_fn = args.input_file + ".clusters.hist.pdf"
        with PdfPages(cl_hist_fn) as pdf_pages:
            for bin_w in [10, 1]:
                try:
                    hist = sns.displot(no_singletons, x="cluster_size", binwidth = bin_w)
                except ValueError:
                    print("There was a problem with the no singletons",
                        f"histogram using binwidth: {bin_w}")
                pdf_pages.savefig()

if __name__ == '__main__':
    desc = ''' Script to summarise the cluster sizes from an MCL clustering'''
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('input_file', nargs='?', metavar='INFILE',
        type=str, default='all-tpm-20-k100.mci.I14', help='Input file name')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
