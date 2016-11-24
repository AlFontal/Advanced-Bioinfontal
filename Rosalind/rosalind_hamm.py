#! /usr/bin/env python
"""
Author: Alejandro Fontal
Student Registration Number: 920110242090
Script
"""
from sys import argv

""" FUNCTIONS """

def calc_ham_dist(seq1, seq2):
    """
    Calculates Hamming distance between two strings of same length.
    :param seq1: String 1
    :param seq2: String 2
    :return: Hamming distance, number of mismatched symbols.
    """
    ham_dist = 0
    for idx, nucleotide in enumerate(seq1):
        if nucleotide != seq2[idx]:
            ham_dist += 1

    return ham_dist

""" Executable code """

if __name__ == "__main__":
    seq_file = open(argv[1], "r")
    seqs = seq_file.read().split("\n")
    seq1 = seqs[0]
    seq2 = seqs[1]
    ham_dist = calc_ham_dist(seq1, seq2)
    output_pm = open('output_pm.txt', "w")
    output_pm.write(str(ham_dist))
    print ham_dist