#!/usr/bin/env python3
#
# Copyright (C) 2015, Universidad Simon Bolivar
# Copying: GNU GENERAL PUBLIC LICENSE Version 2
#
# get_top_predictions.py -- This software scans the semEP predictions
# and then saves them in a file in increasing order.
#
# author: Guillermo Palma (gpalma@ldc.usb.ve)

import sys, subprocess, time, operator

if __name__ == "__main__":
    pred = {}
    fd = open(sys.argv[1], "r")
    for line in fd:
        if "Cluster" != line[:7]:
            tok = line.split("\t")
            key = (tok[0], tok[1])
            prob = float(tok[2][:-1])
            if key in pred:
                if prob > pred[key]:
                    pred[key] = prob
            else:
                pred[key] = prob
    fd.close()
    fd = open(sys.argv[2], "w")
    sorted_pred = sorted(pred.items(), key=operator.itemgetter(1), reverse=True)
    for ((d, t), prob) in sorted_pred:
        fd.write(d+"\t"+t+"\t"+str(prob)+"\n")
    fd.close()
    print("File "+sys.argv[2]+" with the predictions in order was created")

