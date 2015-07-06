#!/usr/bin/env python3
#
# Copyright (C) 2015, Universidad Simon Bolivar
# Copying: GNU GENERAL PUBLIC LICENSE Version 2
#
# test_ding.py -- get the semEP predictions
# on the Yamanishi dataset used in the paper presented by Ding et al [1].
#
# [1] H. Ding, I. Takigawa, H. Mamitsuka, and S. Zhu.
# "Similarity-based machine learning methods for predicting drug-target interactions: A brief review".
# Briefings in Bioinformatics, 2013.
#
# author: Guillermo Palma (gpalma@ldc.usb.ve)

import sys, subprocess, time, operator

nFolds = [str(x) for x in range(1,11)]

def load_interactions(filename):
    interactions = set()
    fd = open(filename, "r")
    for line in fd:
        tok = line.split("\t")
        (d, t) = tok[0], tok[1][:-1]
        interactions.add((d, t))
    return interactions

def check_interactions(dataset):
    interactions = []
    for x in nFolds:
        filename =  "../../test/yamanishi/"+dataset+"/RV_folds/RV_Fold_"+x+".txt"
        inter1 = load_interactions(filename)
        filename =  "../../test/yamanishi/"+dataset+"/folds/OBS_Fold_"+x+".txt"
        inter2 = load_interactions(filename)
        assert( set([]) == inter1 & inter2 )
        interactions.append(inter1 | inter2)
    for i in range(10):
        a = interactions[i]
        for j in range(10):
            b = interactions[j]
            if i != j :
                if a != b:
                    print("Error interactions are diferents")
                    sys.exit(1)
    return interactions[0]

def get_interactions(datasets):
    interactions = {}
    for d in datasets:
        interactions[d] = check_interactions(d)
    return interactions

def get_dataset_pred(sim_d, sim_t, d_interactions):
    pred = {}
    for n in nFolds:
        filename = "Fold_"+n+"_Drugs-Fold_"+n+"_Targets-"+sim_d+"-"+sim_t+"-Pred.txt"
        fd = open(filename, "r")
        for line in fd:
            if "Cluster" != line[:7]:
                tok = line.split("\t")
                key = (tok[0], tok[1])
                prob = float(tok[2][:-1])
                if key not in d_interactions:
                    if key in pred:
                        if prob > pred[key]:
                            pred[key] = prob
                    else:
                        pred[key] = prob
    return pred

def get_predictions(datasets, interactions):
    predictions = {}
    for d in datasets:
        (sim_t, sim_d) = datasets[d]
        predictions[d] = get_dataset_pred(sim_t, sim_d, interactions[d])
    return predictions

def print_prediction(dataset, pred):
    fd = open(dataset+"_predictions.txt", "w")
    sorted_pred = sorted(pred.items(), key=operator.itemgetter(1), reverse=True)
    for ((d, t), prob) in sorted_pred:
        fd.write(d+"\t"+t+"\t"+str(prob)+"\n")

def print_all_predictions(predictions):
    for d in predictions:
        print_prediction(d, predictions[d])

def run_semEP(datasets):
    t_run = []
    for d in datasets:
        start = start_time = time.time()
        for n in nFolds:
            (a, b) = datasets[d]
            semEPcmd = "../../semEP -m -p -t "+a+" -r "+b+" ../../test/yamanishi/"+d+"/sim_mat/Fold_"+n+"_Drugs-Fold_"+n+"_Targets.txt ../../test/yamanishi/"+d+"/"+d+"_desc.txt ../../test/yamanishi/"+d+"/folds_drugs/Fold_"+n+"_Drugs.txt ../../test/yamanishi/"+d+"/folds_targets/Fold_"+n+"_Targets.txt"
            r = subprocess.call(semEPcmd, shell=True)
            if r == 1:
                print("Error executing semEP")
                exit(1)
        tt = time.time() - start_time
        t_run.append((d, tt))
    return t_run

if __name__ == "__main__":
    datasets = {"nr":("0.3421","0.1832"), "gpcr":("0.2759","0.1416"), "ic":("0.2619","0.1355"), "e":("0.2333","0.0209")}
    t_run = run_semEP(datasets)
    interactions = get_interactions(datasets)
    predictions = get_predictions(datasets, interactions)
    print_all_predictions(predictions)
    subprocess.call("rm -rf Fold*", shell=True)
    print("\n***********************************\n")
    print("semEP runtime for all datasets:")
    for (dataset, t) in t_run:
        print('{0} {1:.3f} secs'.format(dataset, t))

    
    
