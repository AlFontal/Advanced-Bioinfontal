#! /usr/bin/env python
"""
Author: Alejandro Fontal
Student Registration Number: 920110242090
Script to solve REVC problem in Rosalind
"""

from sys import argv
fileName = argv[1]

def complement_dna(seq):

    cseq = ""
    seq = seq.upper()
    for nt in seq:
        if nt == "A":
            cseq += "T"
        elif nt == "T":
            cseq += "A"
        elif nt == "C":
            cseq += "G"
        elif nt == "G":
            cseq += "C"

    return cseq[::-1]

if __name__ == "__main__":

    with open(fileName) as seq_file:
        seq = seq_file.read()

        print complement_dna(seq)

