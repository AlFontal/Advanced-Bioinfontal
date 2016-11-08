#!/usr/bin/env python


"""
Author: Alejandro Fontal
Student Registration Number: 920110-242-090
"""


from __future__ import division
from sys import argv


def calc_gc(sequence):
    """

    :param sequence: String containing a DNA sequence. Case insensitive
    :return: Percentage of GC content rounded with 2 decimals.
    """
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    gc_perc = round(gc_count/len(sequence) * 100, 2)

    return gc_perc

def parse_genbank(filename):
    """

    :param filename: Filename of a Genbank file.
    :return: List of dictionaries (and list of tuples) containing Accession ID
    as keys and organism, sequence, GC content (list of tuples) and sequence
    length as values, respectively.
    """
    gb_file = open(filename, "r")
    gb = gb_file.readlines()

    seq = {}
    org = {}
    length = {}
    gc = []
    for idx, line in enumerate(gb):
        if line.startswith("ACCESSION"):
            access = line[12:].split(" ")[0].strip()
            seq[access] = ""
            print access

        elif "ORGANISM" in line:
            organism = line[12:].strip()
            org[access] = organism
            print organism

        elif line.startswith("ORIGIN"):
            seq_idx_start = idx + 1

        elif line.startswith("//"):
            seq_idx_end = idx -1
            seq_lines_len = seq_idx_end - seq_idx_start
            sequence = []
            for i in range(seq_lines_len):
                sequence.append(gb[seq_idx_start + i][10:].strip())
            sequence = "".join(sequence)
            sequence = sequence.replace(" ", "").upper()
            seq[access] = sequence

            gc.append((access, calc_gc(sequence)))

            length[access] = len(sequence)

    parsed_list = [org, seq, gc, length]

    return parsed_list

if __name__ == "__main__":

    abc = parse_genbank("argonaut.gb")

    gc = abc[2]

    print sorted(gc, key=lambda tup: tup[1], reverse= True)

    fasta_file = open("seqs.fasta", "w")

    fasta.file.write






