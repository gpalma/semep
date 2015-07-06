#!/usr/bin/python
#
# Copyright (C) 2014, Universidad Simon Bolivar
# Copying: GNU GENERAL PUBLIC LICENSE Version 2
#
# draw_semEP.py -- Visualization engine for semEP output
# 
# author: Guillermo Palma (gpalma@ldc.usb.ve)

import sys
import random
import os
from pygraphviz import *
import argparse

NODES_PER_CLUSTER = 50
DEFAULT_NODE_SEP = 0.25
NO_PARTITION = -1
OUTPUT_FORMAT = "png"
SEPARATION = 1.5
EDGE_LEN = 3

colors = ["antiquewhite","antiquewhite1","antiquewhite2","antiquewhite3","antiquewhite4","aquamarine","aquamarine1","aquamarine2","aquamarine3","aquamarine4","azure1","azure2","azure3 azure4","beige","bisque","bisque1","bisque2","bisque3","bisque4","black","blue","blue1","blue2","blue3","blue4","blueviolet","brown","brown1","brown2","brown3","brown4","burlywood","burlywood1","burlywood2","burlywood3","burlywood4","cadetblue","cadetblue1","cadetblue2","cadetblue3","cadetblue4","chartreuse","chartreuse1","chartreuse2","chartreuse3","chartreuse4","chocolate","chocolate1","chocolate2","chocolate3","chocolate4","coral","coral1","coral2","coral3","coral4","cornflowerblue","cornsilk","cornsilk1","cornsilk2","cornsilk3","cornsilk4","crimson","cyan","cyan1","cyan2","cyan3","cyan4","darkgoldenrod","darkgoldenrod1","darkgoldenrod2","darkgoldenrod3","darkgoldenrod4","darkgreen","darkkhaki","darkolivegreen","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4","darkorange","darkorange1","darkorange2","darkorange3","darkorange4","darkorchid","darkorchid1","darkorchid2","darkorchid3","darkorchid4","darksalmon","darkseagreen","darkseagreen1","darkseagreen2","darkseagreen3","darkseagreen4","darkslateblue","darkslategray","darkslategray1","darkslategray2","darkslategray3","darkslategray4","darkslategrey","darkturquoise","darkviolet","deeppink","deeppink1","deeppink2","deeppink3","deeppink4","deepskyblue","deepskyblue1","deepskyblue2","deepskyblue3","deepskyblue4","dimgray","dimgrey","dodgerblue","dodgerblue1","dodgerblue2","dodgerblue3","dodgerblue4","firebrick","firebrick1","firebrick2","firebrick3","firebrick4","forestgreen","gold","gold1","gold2","gold3","gold4","goldenrod","goldenrod1","goldenrod2","goldenrod3","goldenrod4","gray","gray0","gray39","gray53","gray78","gray88","green","green1","green2","green3","green4","greenyellow","grey","hotpink","hotpink1","hotpink2","hotpink3","hotpink4","indianred","indianred1","indianred2","indianred3","indianred4","indigo","ivory3","ivory4","khaki","khaki1","khaki2","khaki3","khaki4","lavenderblush2","lavenderblush3","lavenderblush4","lawngreen","lemonchiffon","lemonchiffon1","lemonchiffon2","lemonchiffon3","lemonchiffon4","lightblue","lightblue1","lightblue2","lightblue3","lightblue4","lightcoral","lightcyan","lightcyan1","lightcyan2","lightcyan3","lightcyan4","lightgoldenrod","lightgoldenrod1","lightgoldenrod2","lightgoldenrod3","lightgoldenrod4","lightgoldenrodyellow","lightgray","lightgrey","lightpink","lightpink1","lightpink2","lightpink3","lightpink4","lightsalmon","lightsalmon1","lightsalmon2","lightsalmon3","lightsalmon4","lightseagreen","lightskyblue","lightskyblue1","lightskyblue2","lightskyblue3","lightskyblue4","lightslateblue","lightslategray","lightslategrey","lightsteelblue","lightsteelblue1","lightsteelblue2","lightsteelblue3","lightsteelblue4","lightyellow3","lightyellow4","limegreen","magenta","magenta1","magenta2","magenta3","magenta4","maroon","maroon1","maroon2","maroon3","maroon4","mediumaquamarine","mediumblue","mediumorchid","mediumorchid1","mediumorchid2","mediumorchid3","mediumorchid4","mediumpurple","mediumpurple1","mediumpurple2","mediumpurple3","mediumpurple4","mediumseagreen","mediumslateblue","mediumspringgreen","mediumturquoise","mediumvioletred","midnightblue","mistyrose","mistyrose1","mistyrose2","mistyrose3","mistyrose4","moccasin","navajowhite","navajowhite1","navajowhite2","navajowhite3","navajowhite4","navy","navyblue","olivedrab","olivedrab1","olivedrab2","olivedrab3","olivedrab4","orange","orange1","orange2","orange3","orange4","orangered","orangered1","orangered2","orangered3","orangered4","orchid","orchid1","orchid2","orchid3","orchid4","palegoldenrod","palegreen","palegreen1","palegreen2","palegreen3","palegreen4","paleturquoise","paleturquoise1","paleturquoise2","paleturquoise3","paleturquoise4","palevioletred","palevioletred1","palevioletred2","palevioletred3","palevioletred4","papayawhip","peachpuff","peachpuff1","peachpuff2","peachpuff3","peachpuff4","peru","pink","pink1","pink2","pink3","pink4","plum","plum1","plum2","plum3","plum4","powderblue","purple","purple1","purple2","purple3","purple4","red","red1","red2","red3","red4","rosybrown","rosybrown1","rosybrown2","rosybrown3","rosybrown4","royalblue","royalblue1","royalblue2","royalblue3","royalblue4","saddlebrown","salmon","salmon1","salmon","salmon1","salmon2","salmon3","salmon4","sandybrown","seagreen","seagreen1","seagreen2","seagreen3","seagreen4","seashell3","seashell4","sienna","sienna1","sienna2","sienna3","sienna4","skyblue","skyblue1","skyblue2","skyblue3","skyblue4","slateblue","slateblue1","slateblue2","slateblue3","slateblue4","slategray","slategray1","slategray2","slategray3","slategray4","slategrey","springgreen","springgreen1","springgreen2","springgreen3","springgreen4","steelblue","steelblue1","steelblue2","steelblue3","steelblue4","tan","tan1","tan2","tan3","tan4","thistle","thistle1","thistle2","thistle3","thistle4","tomato","tomato1","tomato2","tomato3","tomato4","turquoise","turquoise1","turquoise2","turquoise3","turquoise4","violet","violetred","violetred1","violetred2","violetred3","violetred4","wheat","wheat1","wheat2","wheat3","wheat4","yellow","yellow1","yellow2","yellow3","yellow4","yellowgreen"]

def load_input_edges(filename):
    edges = set()
    fd = open(filename, "r")
    for line in fd:
        line = line.rstrip("\n")
        tok = line.split("\t")
        if len(tok) == 3 :
            edges.add((tok[0]+" ", tok[1]))
    return edges

def get_name(is_cluster, file_path):
    file_name = os.path.basename(file_path)
    tok = file_name.split("-")
    tname = ""
    if is_cluster:
        n = -4
    else:
        n = -3
    for x in tok[:n]:
            tname += x+"-"
    tname +=tok[n]

    if is_cluster:
        cluster = int(tok[n+1])
    else:
        cluster = NO_PARTITION

    return tname, cluster

def get_nodes_sep(n):
    if n > NODES_PER_CLUSTER :
        return (n/NODES_PER_CLUSTER)*DEFAULT_NODE_SEP
    else:
        return DEFAULT_NODE_SEP

def draw_singleton(root, filename, nodes_singleton, out_format):
    tok = filename.split("-")
    tname = ""
    for x in tok[:-4]:
        tname += x+"-"
    tname +=tok[-4]
    glabel = "\nSingleton nodes: "+tname
    
    Mgraph = AGraph()
    
    # setting graph attributes
    Mgraph.graph_attr["label"] = glabel
    Mgraph.graph_attr["rankdir"] = "LR"
    Mgraph.graph_attr["splines"] = "line"
    Mgraph.graph_attr["overlap"] = "prism"
    Mgraph.graph_attr["bgcolor"] = "ghostwhite"

    # setting node attributes
    Mgraph.node_attr.update(style="filled")
    Mgraph.node_attr.update(fillcolor="white")
    Mgraph.node_attr.update(fontsize="15")

    # Adding singleton nodes
    for node in nodes_singleton:
        Mgraph.add_node(node)
        
    # Draw Partition
    output_file = filename[:-4]+"."+out_format
    print("Partition File: "+output_file)
    out_filepath = os.path.join(root, output_file)
    Mgraph.draw(out_filepath, out_format, "dot")

def get_evidence(node1, node2, evidence_map):
    node1new = node1[:-1]
    evidence = None
    if (node1new, node2) in evidence_map:
        evidence = evidence_map[(node1new, node2)]
    elif (node2, node1new) in evidence_map:
        evidence = evidence_map[(node1new, node2)]
    else:
        print("Error, the pair "+node2+" "+node1new+" in not found in the GO\n")
        sys.exit(1)

    return evidence

def draw_partition(root, filename, nodes1, nodes2, edges, out_format, bp_pred, clean_pred, evidence_map):
    tname, cluster = get_name(True, filename)
    glabel = "\nBipartite Graph: "+tname+"\nPartition: "+str(cluster)
    Mgraph = AGraph()

    # setting graph attributes
    Mgraph.graph_attr["label"] = glabel
    Mgraph.graph_attr["rankdir"] = "LR"
    Mgraph.graph_attr["splines"] = "line"
    Mgraph.graph_attr["overlap"] = "prism"
    Mgraph.graph_attr["bgcolor"] = "ghostwhite"
    ns = get_nodes_sep(len(nodes1)+len(nodes2))
    if evidence_map == None:
        Mgraph.graph_attr["nodesep"] = str(ns)
    else:
        Mgraph.graph_attr["nodesep"] = str(ns*SEPARATION)
    Mgraph.graph_attr["sep"] ="+25,25"

    # setting node attributes
    Mgraph.node_attr.update(style="filled")
    Mgraph.node_attr.update(fillcolor="white")
    Mgraph.node_attr.update(fontsize="15")

    # setting edge attributes
    pos = cluster % len(colors)
    Mgraph.edge_attr.update(color=colors[pos])
    Mgraph.edge_attr.update(style="solid")
    
    # get subgraphs
    gn1 = Mgraph.add_subgraph(nodes1, "subgraph_nodes1")
    gn2 = Mgraph.add_subgraph(nodes2, "subgraph_nodes2")
  
    # adding edges
    for (f, t) in edges:
        if evidence_map == None:
            e = Mgraph.add_edge(f, t)
        else:
            l = get_evidence(f, t, evidence_map)
            size = (len(nodes1)+len(nodes2))/EDGE_LEN
            e = Mgraph.add_edge(f, t, label=l, minlen=str(size))

    # adding predictions
    if len(bp_pred) != 0:  
        Mgraph.edge_attr.update(color="red")
        Mgraph.edge_attr.update(style="dashed")
        for (f, t) in bp_pred:
            e = Mgraph.add_edge(f, t)
   
    if len(clean_pred) != 0:  
        Mgraph.edge_attr.update(color="green")
        Mgraph.edge_attr.update(style="dashed")
        for (f, t) in clean_pred:
            e = Mgraph.add_edge(f, t)
        
    # Draw Partition
    output_file = filename[:-4]+"."+out_format
    print("Partition File: "+output_file)
    out_filepath = os.path.join(root, output_file)
    Mgraph.draw(out_filepath, out_format, "dot")

def load_singleton(file_singleton):
    nodes_singleton = set()
    fd = open(file_singleton, "r")
    for line in fd:
        line = line.rstrip("\n")
        nodes_singleton.add(line)
    return nodes_singleton

def load_cluster(file_cluster):
    fd = open(file_cluster, "r")
    nodes1 = set()
    nodes2 = set()
    edges = set()
    for line in fd:
        line = line.rstrip("\n")
        tok = line.split("\t")
        new_node = tok[0]+" "
        nodes1.add(new_node)
        nodes2.add(tok[1])
        edges.add((new_node, tok[1]))
    return nodes1, nodes2, edges

def get_filepaths(directory):
    file_paths = []
    for root, directories, files in os.walk(directory):
        for filename in files:
            extension = os.path.splitext(filename)[1]
            if extension == ".txt":
                file_paths.append((root, filename)) 
    return file_paths

def get_directory(filename):
    return filename[:-4]+"-Subgr"

def get_predictions(bp_edges, nodes1, nodes2, edges):
    all_edges_partition = set()
    
    for n1 in nodes1:
        for n2 in nodes2:
            all_edges_partition.add((n1, n2))

    all_pred = all_edges_partition - edges
    clean_pred = set()
    bp_pred = set()
    for p in all_pred:
        if p in bp_edges:
            bp_pred.add(p)
        else :
            clean_pred.add(p)
    return bp_pred, clean_pred
    
def draw_and_load_partitions(file_paths, output_format, bp_edges, predictions, evidence_map):
    partitions = []
    all_nodes1 = []
    all_nodes2 = []
    bp_pred = set()
    clean_pred = set()
    nodes_singleton = set()
    for root, filename in file_paths:
        filepath = os.path.join(root, filename)
        if "-singleton.txt" == filepath[-14:] :
            nodes_singleton = load_singleton(filepath)
            draw_singleton(root, filename, nodes_singleton, output_format)
        else:
            nodes1, nodes2, edges = load_cluster(filepath)
            if predictions:
                bp_pred, clean_pred = get_predictions(bp_edges, nodes1, nodes2, edges)
            draw_partition(root, filename, nodes1, nodes2, edges, output_format, bp_pred, clean_pred, evidence_map)
            all_nodes1.extend(nodes1)
            all_nodes2.extend(nodes2)
            partitions.append((filename, edges))
    return partitions, all_nodes1, all_nodes2, nodes_singleton

def draw_principal_graph(inputfile, partitions, nodes1, nodes2, out_format, nodes_singleton):
    tname, cluster = get_name(False, inputfile)
    glabel = "\nBipartite Graph: "+tname
    Mgraph = AGraph()
    
    # setting graph attributes
    Mgraph.graph_attr["label"] = glabel
    Mgraph.graph_attr["rankdir"] = "LR"
    Mgraph.graph_attr["splines"] = "line"
    Mgraph.graph_attr["overlap"] = "prism"
    Mgraph.graph_attr["bgcolor"] = "ghostwhite"
    
    # setting node attributes
    Mgraph.node_attr.update(style="filled")
    Mgraph.node_attr.update(fillcolor="white")
    Mgraph.node_attr.update(fontsize="15")
    
    partition_color = []
    all_edges = []
    for (filename, edges) in partitions:
        tname, cluster = get_name(True, filename)
        pos = cluster % len(colors)
        partition_color.append((colors[pos], pos))
        # setting edge attributes
        for (f, t) in edges:
            ecolor = colors[pos]
            all_edges.append(((f,t),ecolor))  
                       
    # Adding partition information
    cont = 0
    partition_color = sorted(partition_color, key=lambda e: e[1], reverse=True)
    for (ecolor, partition) in partition_color:
        a = cont 
        cont += 1
        b = cont
        cont += 1
        Mgraph.add_node(a, height="0.01", style="invisible", label="")
        Mgraph.add_node(b, height="0.01", style="invisible", label="")
        elabel = "Partition: "+str(partition)
        Mgraph.add_edge(a, b, label=elabel, color=ecolor, minlen="1")

    # Adding singleton nodes
    for node in nodes_singleton:
        Mgraph.add_node(node)

    # Adding partition information
    all_edges.reverse()
    for (f,t), ecolor in all_edges:
        Mgraph.add_edge(f, t, color=ecolor)
        
    # Draw Principal graph
    output_file = inputfile[:-4]+"."+out_format
    print("Output File: "+output_file)
    Mgraph.draw(output_file, out_format, "dot")

def load_go(filename):
    gene_term_map = {}
    fd = open(filename, "r")
    for line in fd:
        tok = line.split("\t")
        gene_term_map[(tok[2], tok[5])] = tok[9]  
        gene_term_map[(tok[0], tok[5])] = tok[9]  
        gene_term_map[(tok[2], tok[4])] = tok[9]  
        gene_term_map[(tok[0], tok[4])] = tok[9]  
    return gene_term_map

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("input_filename", help="File name of the bipartite graph that is the semEP output", action="store")
    parser.add_argument("-o", help="Output format",  action="store",  dest="output_format",  default=OUTPUT_FORMAT)
    parser.add_argument("-p", help="Print preditcted edges",  action="store_true",  dest="prediction",  default=False)
    parser.add_argument("-g", help="Gene Onlology file",  action="store",  dest="go",  default=None)
    args = parser.parse_args()
    
    evidence_map = None
    if args.go != None:
        evidence_map = load_go(args.go)
    random.shuffle(colors)
    bipartite_edges = load_input_edges(args.input_filename)
    directory = get_directory(args.input_filename)
    file_paths = get_filepaths(directory)
    partitions, nodes1, nodes2, nodes_singleton = draw_and_load_partitions(file_paths, args.output_format, bipartite_edges, args.prediction, evidence_map)
    draw_principal_graph(args.input_filename, partitions, nodes1, nodes2, args.output_format, nodes_singleton)

if __name__ == "__main__":
    sys.exit(main())
