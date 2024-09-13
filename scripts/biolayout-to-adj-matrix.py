#! /usr/bin/env python

desc = ''' Script to convert biolayout layout file to an adjacency matrix '''

import argparse
import sys
import re
import numpy as np
import polars as pl
import seaborn as sns

edge = re.compile(r'"(.+)"\s+"(.+)"\s+([0-9\.]+)')

def main(args):
    ''' Main body of code '''

    edge_count_for = {}
    total_edges = 0
    # read in input
    for line in args.layout_file:
        line = line.rstrip('\n')
        if line[0:2] == "//":
            continue
        else:
            match = edge.fullmatch(line)
            if match:
                if args.debug:
                    print(match.group(1, 2, 3))
                total_edges += 1
                for gene in match.group(1, 2):
                    if gene not in edge_count_for:
                        edge_count_for[gene] = 1
                    else:
                        edge_count_for[gene] += 1
            else:
                print("Doesn't match!", line)
    if args.debug:
        print(edge_count_for)

    # print statistics on node degree
    counts = [x for x in edge_count_for.values()]
    print(f"Total Edges = {total_edges}", file = args.summary_file)
    print(f"Mean Edge Degree = {np.mean(counts):.3f}", file = args.summary_file)
    print('Quartiles', file = args.summary_file)
    print(f"Min = {np.quantile(counts, 0):.1f}", file = args.summary_file)
    print(f"Q1 = {np.quantile(counts, 0.25):.1f}", file = args.summary_file)
    print(f"Median = {np.quantile(counts, 0.5):.1f}", file = args.summary_file)
    print(f"Q3 = {np.quantile(counts, 0.75):.1f}", file = args.summary_file)
    print(f"Max = {np.quantile(counts, 1):.1f}", file = args.summary_file)
    print(f"Number of Nodes = {len(edge_count_for.keys())}", file = args.summary_file)

    if args.plot_filebase:
        df = pl.DataFrame({
            "gene": edge_count_for.keys(),
            "num_edges": edge_count_for.values() 
        })
        dens_plot = sns.displot(df, x="num_edges", kind = "kde")
        dens_plot_fn = args.plot_filebase + ".num_edges.dens.png"
        dens_plot.savefig(dens_plot_fn)
        hist_plot = sns.displot(df, x="num_edges", kde = True, binwidth = 50)
        hist_plot_fn = args.plot_filebase + ".num_edges.hist.png"
        hist_plot.savefig(hist_plot_fn)

    if args.output_file:
        for gene in sorted(edge_count_for, key = edge_count_for.get, reverse = True):
            print(f"{gene}\t{edge_count_for[gene]}", file = args.output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('layout_file', nargs='?', metavar='LAYOUT',
        type=argparse.FileType('r'), default='/Users/rjw26/work/apocrita/data/scratch/bty114/zmp-network/zfs-dev-stages-biolayout/publication_change_with_stage_r-0.94.layout.test', 
        help='Biolayout input layout file')
    parser.add_argument('summary_file', nargs='?', metavar='SUMMARY',
        type=argparse.FileType('w'), default=sys.stdout, 
        help='Name of summary stats file')
    parser.add_argument('output_file', nargs='?', metavar='OUTPUT',
        type=argparse.FileType('w'), help='Name of output file')
    parser.add_argument('--plot_filebase', metavar='PLOTS',
        type=str, help='Base name for plot files')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
