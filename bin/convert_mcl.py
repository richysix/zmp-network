#!/usr/bin/env python

desc = ''' Script to parse an mcl graph to nodes and edges file\n

The gene annotation file should contain the following columns\n

node_index - index for node corresponding to current network\n
id - index covering all genes\n
gene_id - Ensembl gene ID\n
gene_name - Gene name\n
description - Description\n
'''
import argparse
import re
import sys

def read_cluster_file(args):
    ''' Opens the cluster file, parses it and returns 2 objects:
    cluster_idx_for, dict: holds the cluster id for each node id
    nodes_for_cluster, dict: holds a list of node ids for each cluster id
    '''
    # colours for nodes
    colours = ["#0073B3", "#CC6600", "#009980", "#F2E640", "#59B3E6", "#CC99B3",
                "#FFFFFF", "#E69900", "#DDDDDD", "#999999"]
    # load cluster file
    nodes = {}
    cluster_size = {}
    nodelist = []
    begin = False
    with open(args.cluster_file, 'r') as f:
        for line in f:
            info = line.lstrip().rstrip('\n').split()
            if not begin:
                if info[0] == "begin":
                    begin = True
                    continue
            else:
                if info[-1] == "$":
                    # remove $
                    info.pop()
                    complete = True
                else:
                    complete = False
                nodelist.extend(info)
                if complete:
                    if args.debug > 1:
                        print("Nodes:", nodelist)
                    cluster_idx = nodelist.pop(0)
                    keep = True
                    if args.cluster_subset:
                        if cluster_idx not in args.cluster_subset:
                            keep = False
                    if len(nodelist) < args.min_cluster_size:
                        keep = False
                    if keep:
                        cluster_size[cluster_idx] = len(nodelist)
                        for node_idx in nodelist:
                            nodes[node_idx] = { 
                                "node_idx": node_idx,
                                "cluster_idx": cluster_idx,
                                "node_colour": colours[int(cluster_idx) % len(colours)]
                            }
                    # reset node list
                    nodelist = []
    return nodes, cluster_size

def read_tab_file(args, nodes):
    ''' Opens the tab file and returns a dict of Ensembl ids for each node id
    '''
    # open tab file
    node_idx_for = {}
    with open(args.graph_tab_file, 'r') as f:
        for line in f:
            node_idx, gene_id = line.rstrip('\n').split('\t') # columns are node_index, gene_id
            if node_idx in nodes:
                nodes[node_idx]['gene_id'] = gene_id
                node_idx_for[gene_id] = node_idx
    return nodes, node_idx_for

def read_annotation_file(args, nodes, node_idx_for):
    ''' Opens the annotation file and returns a dict.
    keys: Ensembl gene ids
    values: a dict containing the annotation information (keys: idx, name, description)
    '''
    with open(args.annotation_file, 'r') as f:
        for idx, line in enumerate(f):
            info = line.rstrip('\n').split('\t') # columns are gene_id, chr, start, end, strand, biotype, name, description
            if info[0] in node_idx_for:
                node_idx = node_idx_for[info[0]]
                nodes[node_idx]['idx'] = str(idx)
                nodes[node_idx]['name'] = info[6]
                nodes[node_idx]['description'] = info[7]
    return nodes

def output_node(args, node, cluster_size, nodes_fh, graphml_fh):
    ''' Open nodes file and output node info as tsv
    '''
    node_idx = node['node_idx']
    cluster_idx = node['cluster_idx']
    node['singleton'] = True if cluster_size[cluster_idx] == 1 else False
    print(node_idx, node["idx"], node["gene_id"], node["name"], 
          node["description"], cluster_idx, str(node['singleton']), 
          sep = "\t", file = nodes_fh)
    if graphml_fh is not None:
        graphml_output_node(args, node, graphml_fh)

def graphml_output_node(args, node, graphml_fh):
    print(f'    <node id="n{node['node_idx']}">', file = graphml_fh)
    print(f'      <data key="d1">{node['idx']}</data>', file = graphml_fh)
    print(f'      <data key="d2">{node['gene_id']}</data>', 
        file = graphml_fh)
    print(f'      <data key="d3">{node['name']}</data>', file = graphml_fh)
    print(f'      <data key="d4">{node['description']}</data>', 
        file = graphml_fh)
    print(f'      <data key="d5">{node['cluster_idx']}</data>', 
        file = graphml_fh)
    print(f'      <data key="d6">{node['singleton']}</data>', 
        file = graphml_fh)
    print(f'      <data key="d7">{node['node_colour']}</data>', 
        file = graphml_fh)
    print(f'    </node>', file = graphml_fh)

def print_graphml_header(args, out_fh):
    graph_id = args.graph_id
    graph_header = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"  
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns 
     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="d0" for="edge" attr.name="weight" attr.type="double"/>
  <key id="d1" for="node" attr.name="node_id" attr.type="int"/>
  <key id="d2" for="node" attr.name="gene_id" attr.type="string"/>
  <key id="d3" for="node" attr.name="name" attr.type="string"/>
  <key id="d4" for="node" attr.name="description" attr.type="string"/>
  <key id="d5" for="node" attr.name="cluster_idx" attr.type="int"/>
  <key id="d6" for="node" attr.name="singleton" attr.type="boolean"/>
  <key id="d7" for="node" attr.name="colour" attr.type="string"/>"""
    print(graph_header, file = out_fh)
    print(f'<graph id="{graph_id}" edgedefault="undirected">', file = out_fh)

def graphml_output_edge(args, edge_idx, source_node_idx, target_node_idx,
                        weight, graphml_fh):
    print(f'    <edge id="e{edge_idx}" source="n{source_node_idx}" target="n{target_node_idx}">',
          file = graphml_fh)
    print(f'      <data key="d0">{weight}</data>', file = graphml_fh)
    print(f'    </edge>', file = graphml_fh)

def read_mci_and_output_edges(args, nodes, graphml_fh):
    # open file and read in input
    begin = False
    edgelist = []
    edge_idx = 1
    edges = set()
    with open(args.edges_file, 'w') as out_fh:
        # add header
        print('edge_idx', 'source', 'target', 'weight', sep="\t", file=out_fh)
        with open(args.mci_file, 'r') as f:
            for line in f:
                info = line.lstrip().rstrip('\n').split()
                if not begin:
                    if info[0] == "begin":
                        begin = True
                        continue
                else:
                    if info[-1] == "$":
                        # remove $
                        info.pop()
                        complete = True
                    else:
                        complete = False
                    edgelist.extend(info)
                    if complete:
                        # go through list and print edges
                        source_node_idx = edgelist.pop(0)
                        source_node_idx = re.sub(':[0-9]+', '', source_node_idx)
                        for edge in edgelist:
                            (target_node_idx, weight) = edge.split(":")
                            # remove self edges
                            if source_node_idx == target_node_idx:
                                continue
                            # check that both nodes exist
                            if source_node_idx in nodes and target_node_idx in nodes:
                                # check whether this edge already exists in opposite direction
                                edge_string = target_node_idx + '-' + source_node_idx
                                if edge_string in edges:
                                    continue
                                else:
                                    edges.add(source_node_idx + '-' + target_node_idx)
                                if args.edge_offset:
                                    weight = float(weight) + args.edge_offset
                                    weight = f"{weight:.3f}"
                                print(edge_idx, source_node_idx, target_node_idx,
                                    weight, sep = "\t", file = out_fh)
                                # output to graphml
                                if graphml_fh is not None:
                                    graphml_output_edge(args, edge_idx, source_node_idx, 
                                        target_node_idx, weight, graphml_fh)
                                edge_idx += 1
                        # clear edge list
                        edgelist = []

def main(args):
    ''' Main body of code '''
    # load cluster file first
    nodes, cluster_size = read_cluster_file(args)

    # load tab file
    nodes, node_idx_for = read_tab_file(args, nodes)

    # load annotation file
    nodes = read_annotation_file(args, nodes, node_idx_for)

    # open output files
    nodes_fh = open(args.nodes_file, 'w')
    print('node_idx', 'id', 'gene_id', 'gene_name', 'description',
          'cluster_id', 'singleton', sep = "\t", file = nodes_fh)
    # open graphML file
    graphml_fh = None
    if args.graphml_file is not None:
        # open file
        graphml_fh = open(args.graphml_file, 'w')
        print_graphml_header(args, graphml_fh)

    # write out nodes file
    for node in nodes:
        if args.debug:
            print(nodes[node])
        output_node(args, nodes[node], cluster_size, nodes_fh, graphml_fh)

    read_mci_and_output_edges(args, nodes, graphml_fh)

    nodes_fh.close()
    if args.graphml_file is not None:
        print("  </graph>\n</graphml>\n", file = graphml_fh)
        graphml_fh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('mci_file', metavar='MCI_FILE',
        type=str, default='graph.mci', help='mcl graph file')
    parser.add_argument('cluster_file', metavar='CLUSTER_FILE',
        type=str, default='graph.mci.I14', help='cluster file from mcl showing which nodes belong to each cluster')
    parser.add_argument('graph_tab_file', metavar='GRAPH_TAB_FILE',
        type=str, default='graph.tab', help='graph tab file containing mcl node index and gene id')
    parser.add_argument('annotation_file', metavar='ANNOTATION_FILE',
        type=str, default='annotation.txt', help='annotation file containing gene id and name')
    parser.add_argument('--graphml_file', metavar='GRAPHML FILE',
        type=str, default=None, help='Name of the graphml output file if required')
    parser.add_argument('--nodes_file', metavar='NODES FILE',
        type=str, default='graph.nodes.tsv', help='Name of the output file for the nodes')
    parser.add_argument('--edges_file', metavar='EDGES FILE',
        type=str, default='graph.edges.tsv', help='Name of the output file for the edges')
    parser.add_argument('--edge_offset', metavar='EDGE OFFSET',
        type=float, default=None, help='An amount to add to edge weights')
    parser.add_argument('--graph_id', metavar='GRAPH ID',
        type=str, default="G", help='Id for the graphml file')
    parser.add_argument('--min_cluster_size', metavar='MINIMUM CLUSTER SIZE',
        type=int, default=2, help='The minimum size for clusters to keep')
    parser.add_argument('--cluster_subset', metavar='CLUSTER SUBSET',
        type=str, default=[], action = "append",
        help='Cluster to subset output to. Can be specified multiple times')
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
# This software is Copyright (c) 2024. Queen Mary, University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
