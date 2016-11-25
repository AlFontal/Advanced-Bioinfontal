#! /usr/bin/env python
"""
Author: Alejandro Fontal
Student Registration Number: 920110242090
Script to solve HAMM problem in Rosalind
"""

from sys import argv


def calc_ham_dist(seq1, seq2):
    """
    Calculates Hamming distance between two strings of same length.
    :param seq1: String 1
    :param seq2: String 2
    :return: Hamming distance, number of mismatched symbols.
    """
    ham_dist = 0
    seqs = zip(seq1, seq2)
    for nt in seqs:
        if nt[0] != nt[1]:
            ham_dist += 1

    return ham_dist


if __name__ == "__main__":

    with open(argv[1]) as seq_file:
        seq1, seq2 = seq_file.read().split("\n")

    ham_dist = calc_ham_dist(seq1, seq2)
    with open("output_pm.txt", "w") as output_pm:
        output_pm.write(str(ham_dist))

    print ham_dist
