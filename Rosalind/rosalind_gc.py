#! /usr/bin/env python
"""
Author: Alejandro Fontal
Student Registration Number: 920110242090
Script to solve GC problem in Rosalind
"""


from __future__ import division
from sys import argv
fileName = argv[1]


def parse_fasta(fasta_fn):
    """
    Takes a fasta file and parses it returning a dictionary containing
    labels as keys and sequences as values.
    :param fasta_fn: Filename of the FASTA file to parse
    :return: A dictionary containing labels as keys and sequences as values.
    """
    with open(fasta_fn) as fasta_file:
        fasta_list = fasta_file.read().splitlines()
        parsed_seqs = {}
        for line in fasta_list:
            if line.startswith(">"):
                label = line[1:]
                parsed_seqs[label] = ""
            else:
                parsed_seqs[label] += line
    return parsed_seqs


def calc_max_GC(seqs_dict):
    """
    Takes a list of sequences and gives the GC% of the one with the most GC%

    :param seqs: List of DNA sequences in string form
    :param names: List of DNA sequence names.
    :return: Two values. The maximum percentage of GC content of the sequences
    and the name of the sequence with the maximum percentage.
    """

    max_gc = 0

    for seq in seqs_dict.items():
        label, seq = seq
        gc = seq.count("G") + seq.count("C")
        gc_perc = (gc/len(seq)) * 100
        if gc_perc > max_gc:
            max_gc = gc_perc
            max_label = label

    return max_label, max_gc


if __name__ == "__main__":

    seqs = parse_fasta(fileName)

    max_label, max_gc = calc_max_GC(seqs)

    with open('output.txt', 'w') as f:
        out_str = "{}\n{}".format(max_label, max_gc)
        f.write(out_str)

    print out_str