#! /usr/bin/env python

from __future__ import division
from sys import argv
import numpy as np


"""
Author: Alejandro Fontal and Roos Goessen
Script to solve P3
"""


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


def fastq_stats(normal_dict, trimmed_dict):
    """

    :param fastq_dict:
    :return:
    """

    """
    Calculate Stats for the non-trimmed data
    """

    lengths1 = map(len, normal_dict.keys())
    max_length1 = max(lengths1)
    min_length1 = min(lengths1)
    avg_length1 = np.mean(lengths1)

    qual_dict_1 = {}

    for read in normal_dict:
        quality_list = normal_dict[read]
        for idx, value in enumerate(quality_list):
            if idx in qual_dict_1.keys():
                qual_dict_1[idx].append(value)
            else:
                qual_dict_1[idx] = [value]

    for i in range(max_length1):
        qual_dict_1[i] = round(np.mean(qual_dict_1[i]), 2)


    """
    Calculate Stats for the trimmed data
    """


    lengths2 = map(len, trimmed_dict.keys())
    max_length2 = max(lengths2)
    min_length2 = min(lengths2)
    avg_length2 = np.mean(lengths2)

    qual_dict_2 = {}

    for read in trimmed_dict:
        quality_list = trimmed_dict[read]
        for idx, value in enumerate(quality_list):
            if idx in qual_dict_2.keys():
                qual_dict_2[idx].append(value)
            else:
                qual_dict_2[idx] = [value]

    for i in range(max_length2):
        qual_dict_2[i] = round(np.mean(qual_dict_2[i]), 2)

    diff = []
    for i in range(len(qual_dict_1)):
        diff.append(round(qual_dict_2[i]-qual_dict_1[i], 3))




    print "ORIGINAL:    min = {}    max = {}    " \
          "avg = {}".format(max_length1, min_length1, avg_length1)

    print "TRIMMED:    min = {}    max = {}    " \
          "avg = {}".format(max_length2, min_length2, avg_length2)

    print "The average read quality per position is:\n "
    print "{}\t\t{}\t\t{}\t\t{}\n".format("Position".center(15),
                            "Original".center(15), "Trimmed".center(15),
                                    "Difference".center(15))
    for i in range(max_length1):
        print "{}\t\t{}\t\t{}\t\t{}".format(str(i+1).center(15),
                                            str(qual_dict_1[i]).center(15),
                                            str(qual_dict_2[i]).center(15),
                                            str(diff[i]).center(15))

    return qual_dict_1




if __name__ == "__main__":

    #filename = argv[1]
    filename = "tomatosample.fq"
    seqs_dict = parse_fastq_file(filename)


    filename = "trimmed.fq"
    seqs_dict2 = parse_fastq_file(filename)
    trimmed_dic = fastq_stats(seqs_dict, seqs_dict2)




