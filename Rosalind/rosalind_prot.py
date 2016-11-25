#! /usr/bin/env python

from sys import argv

"""
Author: Alejandro Fontal
Student Registration Number: 920110242090
Script to solve REVC problem in Rosalind
"""

def get_codons_dict(codons_table):

    codons_dict = {}
    items = codons_table.split()
    for idx, item in enumerate(items):
        if len(item) == 3:
            codons_dict[item] = items[idx + 1]

    return codons_dict

def rna_to_prot(seq, codons_dict):

    prot = ""
    for i in range(0, len(seq)-2, 3):

        if codons_dict[seq[i:i+3]] == "Stop":
            return prot

        else:
            prot += codons_dict[seq[i:i+3]]



if __name__ == "__main__":

    with open("codons_table.txt") as codons_file:
        codons_table = codons_file.read()

    with open(argv[1]) as seq_file:
        seq = seq_file.read()

    codons_dict = get_codons_dict(codons_table)

    prot = rna_to_prot(seq, codons_dict)
    print prot
