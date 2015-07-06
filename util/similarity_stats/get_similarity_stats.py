#!/usr/bin/env python
#
# Copyright (C) 2014, Universidad Simon Bolivar
# Copying: GNU GENERAL PUBLIC LICENSE Version 2
#
# get_similarity_stats.py -- Simple similarity measure stats
# 
# author: Guillermo Palma (gvpalma@usb.ve)

import sys, os, traceback, numpy

def load_sim(filename):
    filename = open(filename, "r")
    sim = []
    for line in filename:
        tok = line.split("\t")
        if tok[0] != tok[1] :
            sim.append( float(tok[2][:-1] ) )
    return sim

def get_stats(a):
    print("Similarity measure stats:")
    print("Min: {0:.4f}".format(numpy.amin(a)))
    print("Max: {0:.4f}".format(numpy.amax(a)))
    print("Average: {0:.4f}".format(numpy.average(a)))
    print("Median: {0:.4f}".format(numpy.median(a)))
    assert(numpy.percentile(a, 50) ==  numpy.median(a))
    print("****************************")
    print("Percentile\tValue\n")
    percentiles = [x*10 for x in range(1,10)]
    for p in percentiles:
        print("{0}\t{1:.4f}".format(p, numpy.percentile(a, p)))

if __name__ == "__main__":
    sim = load_sim(sys.argv[1])
    get_stats(sim)
 
