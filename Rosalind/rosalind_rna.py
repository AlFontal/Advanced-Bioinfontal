#! /usr/bin/env python
"""
Author: Alejandro Fontal
Student Registration Number: 920110242090
Script to solve RNA problem in Rosalind
"""


from __future__ import division
from sys import argv
fileName = argv[1]

def dna_to_rna(seq):
    """
    Given a string, returns same string with T's replaced by Us

    :param seq: String
    :return: String with T's replaces by U's
    """

    seq = seq.upper()
    seq = seq.replace("T", "U")

    return seq

if __name__ == "__main__":

    with open(fileName) as seq_file:
        seq = seq_file.read()
        print dna_to_rna(seq)
