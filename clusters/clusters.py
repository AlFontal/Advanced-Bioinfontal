#!/usr/local/bin/python

import numpy
import scipy
import scipy.cluster.hierarchy as sch
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import subprocess
"""
Script to practise clustering of biological data
"""
__author__ = 'Alejandro Fontal'
__email__ = "alejandro.fontal@wur.nl"

def parse_txt(txt_fn):
    """
    Parses a tab delimited text file and returns a list with column headers,
    a list with row headers and a dictionary with first element of each row
    as key and the rest of the elements of the row as values
    (will only work if they are floats).

    :param txt_fn: Filename of the txt file to parse
    :return: genes (dict) and labels (list)
    """
    txt_file = open(txt_fn, "r")
    lines = txt_file.readlines()
    samples = lines[0].strip().split("\t")[1:]
    data = lines[1:]
    exps = {}
    genes = []
    for line in data:
        line = line.split("\t")
        gene = line[0]
        exp_vals = map(float, line[1:])
        exps[gene] = exp_vals
        genes.append(gene)
    return samples, genes, exps

if __name__ == "__main__":

    txt_fn = "MayAppleExpressionData.txt"
    samples, genes, exps = parse_txt(txt_fn)

    data_frame = pd.DataFrame(exps, index=samples)
    distances = sch.distance.pdist(data_frame, metric="euclidean")
    clustering = sch.linkage(distances, method='complete')
    tree = sch.dendrogram(clustering,labels=samples)
    plt.show()

    data_frame = data_frame.transpose()

    distances = sch.distance.pdist(data_frame, metric = 'correlation')
    clustering = sch.linkage(distances, method='complete')
    tree = sch.dendrogram(clustering , leaf_font_size=2,
                          color_threshold=4, labels= -genes)
    plt.show()

    clustermap_fn = "testCluster.pdf"
    clustermap = sns.clustermap(data_frame.transpose(),
                                figsize=(20,12),
                                metric= "correlation",
                                method = "complete",
                                row_cluster=False)
    clustermap.savefig(clustermap_fn)
    cmd = "evince {}".format(clustermap_fn)
    subprocess.check_call(cmd, shell=True)



