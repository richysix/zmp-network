#! /usr/bin/env python

''' Script to subset a correlation network to a specific cluster '''

import argparse
import sys
import csv

def main(args):
    ''' Main body of code '''

    # open nodes output file
    with open(args.nodes_outfile, 'w', newline='') as csvfile:
        nodes_out = csv.writer(csvfile, quotechar='"', delimiter=',',
                                quoting=csv.QUOTE_ALL)
 
        # read in input
        nodes = []
        header = True
        lines = csv.reader(args.nodes_file, quotechar='"', delimiter=',',
                        quoting=csv.QUOTE_ALL, skipinitialspace=True)
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

    with open(args.edges_outfile, 'w', newline='') as csvfile:
        edges_out = csv.writer(csvfile, delimiter=',',
                                quoting=csv.QUOTE_MINIMAL)
        # read in edges file
        header = True
        lines = csv.reader(args.edges_file, quotechar='"', delimiter=',',
                        quoting=csv.QUOTE_ALL, skipinitialspace=True)
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
    parser.add_argument('nodes_file', nargs='?', metavar='FILE',
        type=argparse.FileType('r'), default=sys.stdin, help='')
    parser.add_argument('edges_file', nargs='?', metavar='FILE',
        type=argparse.FileType('r'), default=sys.stdin, help='')
    parser.add_argument('--cluster_num',
        type=str, default="1", help='Number of cluster to subset to')
    parser.add_argument('--nodes_outfile',
        type=str, default="nodes.csv", help='Name of file to write the nodes to')
    parser.add_argument('--edges_outfile',
        type=str, default="edges.csv", help='Name of file to write the edges to')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
