#! /usr/bin/env python

desc = ''' Script to create counts for a histogram of correlation values '''

import argparse
import os
import csv
import math
import sys
import re
import gzip

def main(args):
    ''' Main body of code '''
    # work out delimiter based on filename extension
    filename, extension = os.path.splitext(args.cor_mat_file)
    gzipped = False
    if extension == ".gz":
        gzipped = True
        filename, extension = os.path.splitext(filename)
    delim = "," if extension == ".csv" else "\t"
    hist_counts = {}
    for i in range(-10, 10, 1):
        j = str(i + 1)
        i = str(i)
        hist_counts[(i, j)] = 0
    if args.debug > 1:
        print(hist_counts)
    
    # read in input
    if gzipped:
        csv_in = gzip.open(args.cor_mat_file, 'rt', newline='')
    else:
        csv_in = open(args.cor_mat_file, 'r', newline='')

    # read in input
    file_lines = csv.reader(csv_in, quotechar='"', delimiter=delim,
                        quoting=csv.QUOTE_ALL, skipinitialspace=True)
    header = True
    for line in file_lines:
        if header:
            header = False
            continue
        for cor in line[1:]:
            if abs(float(cor)) == 1:
                continue
            lower = math.floor(float(cor)*10)
            upper = math.ceil(float(cor)*10)
            if lower == upper:
                lower = upper - 1
            hist_counts[(str(lower), str(upper))] += 1

    if args.debug > 1:
        print(hist_counts)
    for k in hist_counts.keys():
        start, stop = k
        print(f"{int(start)/10:.1f}\t{int(stop)/10:.1f}\t{hist_counts[k]}",
            file = args.output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cor_mat_file', metavar='COR_MATRIX_FILE',
        type=str, default='all-tpm-orig.mat.csv', 
        help='File with a matrix of cor values')
    parser.add_argument('output_file', nargs='?', metavar='OUTFILE',
        type=argparse.FileType('w'), default=sys.stdout, help='Output file name')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    args = parser.parse_args()
    main(args)

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2024. University of Cambridge.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
