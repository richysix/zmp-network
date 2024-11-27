#! /usr/bin/env python

''' Script to take a DETCT counts file, subset to required samples
and aggregate the count data from regions to genes '''

import argparse
import os
import csv

def main(args):
    ''' Main body of code '''
    nodes_filename, nodes_extension = os.path.splitext(args.nodes_file)
    nodes_delim = "," if nodes_extension == ".csv" else "\t"
    
    nodes = []
    nodes_for = {}
    with open(args.nodes_file, 'r', newline='') as csv_in:
        # read in input
        header = True
        lines = csv.reader(csv_in, quotechar = '"', delimiter = nodes_delim,
                            quoting = csv.QUOTE_ALL, skipinitialspace = True)
        for line in lines:
            if header:
                header = False
                continue
            if line[5] in nodes_for:
                nodes_for[line[5]].append(line)
            else:
                nodes_for[line[5]] = [line]
        
    # open file for all selected nodes
    all_file = args.out_file_base + ".all.tsv"
    with open(all_file, 'w', newline = '') as allfile:
        all_nodes = csv.writer(
            allfile, quotechar='"', delimiter="\t", 
            quoting=csv.QUOTE_MINIMAL
        )
        for cluster_id in nodes_for:
            # only output nodes to cluster file
            # from clusters bigger than the minimum size
            if len(nodes_for[cluster_id]) >= args.min_cluster_size:
                output = True
                out_file = args.out_file_base + ".cluster-" + str(cluster_id) + ".tsv"
                genes_file = open(out_file, 'w')
            else:
                output = False
            for node in nodes_for[cluster_id]:
                all_nodes.writerow(node)
                if output:
                    print(node[2], file = genes_file)
            if output:
                genes_file.close()

if __name__ == '__main__':
    desc = ''' Script to take a nodes files and output an all genes file and
a file for the genes for each cluster'''
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('nodes_file', metavar='INFILE',
        type=str, default='all.csv', help='Input file name')
    parser.add_argument('out_file_base', metavar='OUTFILE',
        type=str, default='all-tpm', help='Base name for the output files')
    parser.add_argument('--min_cluster_size', default=10, type=int,
        help='Minimum size of cluster')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    if params.debug > 1:
        print(params)
    main(params)
