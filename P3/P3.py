#! /usr/bin/env python
from __future__ import division

"""
Author: Alejandro Fontal and Roos Goessen
Script to solve P3
"""

from sys import argv



def parse_fastq_file(filename):
    """

    :param filename: Filename of a fastQ file with Illumina 1.5+ Encoding
    :return: A dictionary with sequences as values and quality scores ranging
    from 0 to 41.
    """
    filename = open(filename, "r")
    file = filename.readlines()

    seqs_dict = {}
    for line in file:
        if line.startswith("@"):
            seq = True
        elif seq == True:
            curr_seq = line.strip()
            seqs_dict[curr_seq] = []
            seq = False
        elif line.startswith("+"):
            qc = True
            continue
        elif qc == True:
            qc_vals = []
            for char in line.strip():
                qual = ord(char) - 64
                qc_vals.append(qual)
                seqs_dict[curr_seq] = qc_vals
            qc = False

    return seqs_dict

def fastq_stats(fastq_dict):
    """

    :param fastq_dict:
    :return:
    """
    lengths = map(len, fastq_dict.keys())
    max_length = max(lengths)
    min_length = min(lengths)
    avg_length = mean(lengths)







if __name__ == "__main__":

    seqs_dict = parse_fastq_file("tomatosample.fq")

    lengths = map(len, seqs_dict.keys())
    print max(lengths)