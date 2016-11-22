#!/usr/local/bin/python

from __future__ import division
import subprocess

"""
Script to solve Exam from Nov'15 in Advanced Bioinformatics
"""
__author__ = 'Alejandro Fontal'
__email__ = "alejandro.fontal@wur.nl"


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


if __name__ == "__main__":

    chrom = parse_fasta("chr3.fa")
    contigs = parse_fasta("velvet_15.fa")
