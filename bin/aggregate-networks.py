#! /usr/bin/env python

''' Script to aggregate mulitple coexpression networks '''

import argparse
import sys
import numpy as np
from scipy.stats import rankdata
import pandas as pd
from random import sample
import gzip

def load_matrix_file(filename):
    net = pd.read_csv(filename, sep="\t", header=0, dtype = {"dummy": "str"})
    idx = net["dummy"]
    net = net.drop(columns = ["dummy"])
    tab_file = filename.replace("mat", "tab")
    tab = pd.read_csv(tab_file, sep="\t", header=None, names=("idx", "GeneID"), index_col = 0)
    net.columns = tab.loc[[int(x) for x in idx], "GeneID"]
    net.index = tab.loc[[int(x) for x in idx], "GeneID"]
    return(net, tab)

def rank_and_standardise(matrix, absolute=False):
    cols = matrix.columns
    rows = matrix.index
    matrix = matrix.to_numpy()
    np.fill_diagonal(matrix, 0)
    if (absolute):
        matrix = np.absolute(matrix)
    matrix = rankdata(matrix, method="min").reshape(matrix.shape)
    matrix = matrix/matrix.max()
    return(
        pd.DataFrame(
            matrix,
            index = rows,
            columns = cols
        )
    )

def aggregate_networks(net1, net2, tab):
    cols = [ x for x in tab["GeneID"] ]
    net2 = pd.DataFrame(net2, index = cols, columns = cols)
    # add net2 subset to correct indices of net1
    net1.loc[cols, cols] = net1.loc[cols, cols] + net2
    return(net1)

def load_annotation(anno_file, args = None):
    idx = 0
    idx_for = {}
    gene_id_for = {}
    with open(anno_file) as anno_fh:
        with open(f"{args.output_base}-network.nodes.tsv", mode='w') as nodes_fh:
            print(f"node_idx\tgene_id", file=nodes_fh)
            col_name_for = {
                0: "GeneID",
                1: "Chr",
                2: "Start",
                3: "End",
                4: "Strand",
                5: "Biotype",
                6: "Name",
                7: "Description"
            }
            for line in anno_fh:
                fields = line.rstrip('\n').split("\t")
                info = {col_name_for[i]: field for i, field in enumerate(fields)}
                if args and args.debug > 1:
                    print(info)
                print(f"{idx}\t{info['GeneID']}", file=nodes_fh)
                idx_for[info['GeneID']] = idx
                gene_id_for[idx] = info['GeneID']
                idx += 1
    return(idx_for, gene_id_for)

def output_network(net, filename):
    # set diagonal to zero
    np.fill_diagonal(net.to_numpy(), 0)
    # remove any rows/columns with all zeros
    to_keep = list(net.sum(0) != 0)
    net = net.loc[to_keep, to_keep]
    net.to_csv(
        path_or_buf = filename,
        sep = "\t",
        header = True
    )

# @profile
def main(args):
    ''' Main body of code '''
    # load gene annotation
    idx_for, gene_id_for = load_annotation(args.annotation, args)
    if args.debug > 1:
        print(idx_for, file=sys.stderr)
        print(gene_id_for, file=sys.stderr)
    num_genes = len(idx_for)
    fofn = [sample(args.input_files, len(args.input_files)) for x in range(args.orderings)]
    print("Files:\n", file=sys.stderr)

    # loop through input files
    for i, files in enumerate(fofn):
        print(files, file=sys.stderr)
        # make base aggregate network
        agg = pd.DataFrame(
            np.zeros(shape=(num_genes, num_genes)),
            index = gene_id_for.values(),
            columns = gene_id_for.values()
        )
        for j, file in enumerate(files):
            # read in matrix file
            net, tab = load_matrix_file(file)
            # rank and standardise
            net = rank_and_standardise(net, absolute=True)
            # add to aggregate
            agg = aggregate_networks(agg, net, tab)
            # rank and standardise
            agg = rank_and_standardise(agg)
            # output aggregated network
            filename = f"{args.output_base}-ordering{i+1}-networks{j+1}.mat.gz"
            output_network(agg, filename)
        filename = f"{args.output_base}-ordering{i+1}-networks{j+2}.mat.gz"
        output_network(net, filename)
        filename = f"{args.output_base}-ordering{i+1}-networks{j+2}.edges.gz"
        # prune network to top 10%
        threshold = np.quantile(agg, 0.9)
        agg[agg < threshold] = 0
        # print network as edges file
        # as most edges have been set to zero
        colnames = list(agg.columns.values)
        with gzip.open(filename, 'wt') as out_fh:
            print(f"edge_idx\tsource\ttarget\tweight", file = out_fh)
            edge_i = 0
            for i in range(len(colnames)):
                row = colnames[i]
                j = 0
                while i > j:
                    col = colnames[j]
                    if agg.loc[row, col] != 0:
                        print(f"{edge_i}\t{idx_for[row]}\t{idx_for[col]}\t{agg.loc[row, col]}", file = out_fh)
                        edge_i += 1
                    j += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_files', nargs='*', metavar='FILE',
        type=str, default=sys.stdin, 
        help='Input matrix files')
    parser.add_argument('--output_base', metavar='OUT FILE BASE',
        type=str, default="agg", 
        help='Base file name for the output networks')
    parser.add_argument('--annotation', metavar='ANNOTATION FILE',
        type=str, default="annotation.txt", 
        help='Gene annotation file')
    parser.add_argument('--orderings', metavar='INT',
        type=int, default=10, 
        help='Number of different file orderings to do')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
