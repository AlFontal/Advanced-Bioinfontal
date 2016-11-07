#! /usr/bin/env python
from __future__ import division

"""
Author: Alejandro Fontal and Roos Goessen
Script to solve P3
"""

from sys import argv



def parse_fastq_file(filename):

    fastq_dict = {}
    line_pos = 0
    filename = filename.split()
    value_line = []
    for line in filename:
        print line
        if line.startswith('@'):
            label = line
            fastq_dict[label] = ""
            line_pos = 1
        elif line_pos == 1:
            for i in line:
                value = ord(i)
                value = (((value-33)/93)*41)
                value = '%.f' % value
                value = str(value)
                value_line.append(value)

            fastq_dict[label] = value_line
            value_line = []

    return fastq_dict


if __name__ == "__main__":

    filename = open("tomatosample.fq", "r")
    file = filename.readlines()

    seqs = {}
    for line in file:
        if line.startswith("@"):
            seq = True
        elif seq == True:
            curr_seq = line.strip()
            seqs[curr_seq] = []
            seq = False
        elif line.startswith("+"):
            qc = True
            continue
        elif qc == True:
            qc_vals = []
            for char in line.strip():
                qual = ord(char) - 64
                qc_vals.append(qual)
                seqs[curr_seq] = qc_vals
            qc = False

    print seqs.items()[0:5]


    #dict = parse_fastq_file(file)

    #print dict['@FCC0U42ACXX:2:1101:11738:4487#ACTACAAG/1']

