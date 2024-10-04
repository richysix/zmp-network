#! /usr/bin/env python

''' Script to subset a correlation network to a specific cluster '''

import argparse
import sys
import csv
import os

def main(args):
    ''' Main body of code '''

    # work out delimiter based on filename extension
    nodes_filename, nodes_extension = os.path.splitext(args.nodes_file)
    nodes_delim = "," if nodes_extension == ".csv" else "\t"
    edges_filename, edges_extension = os.path.splitext(args.edges_file)
    edges_delim = "," if edges_extension == ".csv" else "\t"
    # create output filenames if not specified
    if args.nodes_outfile is None:
        args.nodes_outfile = ".".join([nodes_filename, "cluster" + args.cluster_num, nodes_extension])
    if args.edges_outfile is None:
        args.edges_outfile = ".".join([edges_filename, "cluster" + args.cluster_num, edges_extension])
    if args.debug:
        print(f"Output filenames: nodes = {args.nodes_outfile}, edges = {args.edges_outfile}")

    # open nodes input file
    nodes = []
    with open(args.nodes_file, 'r', newline='') as csv_in:
        # read in input
        header = True
        lines = csv.reader(csv_in, quotechar='"', delimiter=nodes_delim,
                            quoting=csv.QUOTE_ALL, skipinitialspace=True)
        # open nodes output file
        with open(args.nodes_outfile, 'w', newline='') as csvfile:
            nodes_out = csv.writer(csvfile, quotechar='"', delimiter=nodes_delim,
                                    quoting=csv.QUOTE_ALL)
            for line in lines:
                if header:
                    header = False
                    nodes_out.writerow(line)
                    continue
                if args.debug:
                    print(line)
                if line[5] == args.cluster_num:
                    nodes.append(line[0])
                    nodes_out.writerow(line)

    with open(args.edges_file, 'r', newline='') as csv_in:
        # read in edges file
        header = True
        lines = csv.reader(csv_in, quotechar='"', delimiter=edges_delim,
                        quoting=csv.QUOTE_ALL, skipinitialspace=True)
        with open(args.edges_outfile, 'w', newline='') as csvfile:
            edges_out = csv.writer(csvfile, delimiter=edges_delim,
                                    quoting=csv.QUOTE_MINIMAL)
            for line in lines:
                if header:
                    header = False
                    edges_out.writerow(line)
                    continue
                if args.debug:
                    print(line)
                if line[1] in nodes and line[2] in nodes:
                    edges_out.writerow(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('nodes_file', metavar='NODES FILE',
        type=str, default='nodes.csv', 
        help='Name of the nodes file for the whole network')
    parser.add_argument('edges_file', metavar='EDGES FILE',
        type=str, default='edges.csv', 
        help='Name of the edges file for the whole network')
    parser.add_argument('--cluster_num',
        type=str, default="1", help='Number of cluster to subset to')
    parser.add_argument('--nodes_outfile',
        type=str, default=None, help='Name of file to write the nodes to. If not specified, created from input file name and cluster number')
    parser.add_argument('--edges_outfile',
        type=str, default=None, help='Name of file to write the edges to. If not specified, created from input file name and cluster number')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
