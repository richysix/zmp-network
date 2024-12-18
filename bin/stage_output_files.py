#! /usr/bin/env python

''' Script to select files for staging '''

import argparse
import sys
from pathlib import Path

def main(args):
    ''' Main body of code '''

    # read in input
    header = True
    for line in args.input_file:
        items = line.rstrip('\n').split("\t")
        if header:
            col_idx_for = { items[x]: x for x in range(len(items)) }
            header = False
            continue
        if args.expt:
            if items[col_idx_for["expt"]] != args.expt:
                continue
        if args.method:
            if items[col_idx_for["method"]] != args.method:
                continue
        info = { name: items[col_idx_for[name]] for name in ["expt", "threshold", "inflation"] }
        # create new dir
        Path(info['expt']).mkdir(exist_ok = True)
        dir = Path('.')
        # filtered network
        file = Path(f"{info['expt']}-all-tpm-t20-k{info['threshold']}.mcx")
        new_file = dir / info['expt'] / file.name
        file.rename(new_file)

        filename = f"{info['expt']}-all-tpm-t20-k{info['threshold']}.mcx.I{info['inflation']}"
        # get filenames matching pattern
        for file in dir.glob(f"{filename}*"):
            # try moving file it to new dir
            new_file = dir / info['expt'] / file.name
            file.rename(new_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', nargs='?', metavar='FILE',
        type=argparse.FileType('r'), default=sys.stdin, 
        help='Name of file with info on which files to select')
    parser.add_argument('--expt', metavar='FILE', type=str, default=None,
        help='Name of an experiment to limit to')
    parser.add_argument('--method', metavar='FILE', type=str, default=None,
        help='Name of an experiment to limit to')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
