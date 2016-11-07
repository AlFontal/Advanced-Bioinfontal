#! /usr/bin/env python
from __future__ import division

"""
Author: Alejandro Fontal and Roos Goessen
Script to solve P3
"""

from sys import argv
import numpy as np


def parse_fastq_file(filename):
    """

    :param filename: Filename of a fastQ file with Illumina 1.5+ Encoding
    :return: A dictionary with sequences as values and quality scores ranging
    from 0 to 41.
    """
    filename = open(filename, "r")
    file = filename.readlines()

    seqs_dict = {}
    is_seq = False
    is_qc = False
    for line in file:
        if line.startswith("@"):
            is_seq = True

        elif is_seq:
            curr_seq = line.strip()
            seqs_dict[curr_seq] = []
            is_seq = False

        elif line.startswith("+"):
            is_qc = True

        elif is_qc:
            qc_vals = []
            for char in line.strip():
                qual = ord(char) - 64  # Convert Illumina scores to 0-41 scores
                qc_vals.append(qual)
                seqs_dict[curr_seq] = qc_vals
            is_qc = False

    return seqs_dict

def fastq_stats(fastq_dict):
    """

    :param fastq_dict:
    :return:
    """
    lengths = map(len, fastq_dict.keys())
    max_length = max(lengths)
    min_length = min(lengths)
    avg_length = np.mean(lengths)

    qual_dict = {}

    for read in fastq_dict:
        quality_list = fastq_dict[read]
        for idx, value in enumerate(quality_list):
            if idx in qual_dict.keys():
                qual_dict[idx].append(value)
            else:
                qual_dict[idx] = [value]

    for i in range(max_length):
        qual_dict[i] = round(np.mean(qual_dict[i]), 2)

    print "\nThe maximum sequence length is {}\n".format(max_length)
    print "The minimum sequence length is {}\n".format(min_length)
    print "The average sequence length is {}\n".format(avg_length)
    print "The average read quality per position is:\n "
    print "{}  {}\n".format("Position".center(15),
                            "Average Quality".center(15))
    for i in range(max_length):
        print "{}     {}".format(str(i+1).center(15),
                                 str(qual_dict[i]).ljust(15))






if __name__ == "__main__":

    #filename = argv[1]
    filename = "tomatosample.fq"
    seqs_dict = parse_fastq_file(filename)
    fastq_stats(seqs_dict)
