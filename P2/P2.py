#!/usr/bin/env python

from __future__ import division

"""
Author: Alejandro Fontal
Student Registration Number: 920110-242-090
"""


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
    if type(filename) == str:
         gb_file = open(filename, "r")
         gb = gb_file.readlines()

    if type(filename) == list:
        gb = filename

    else:
        print "Input file type not supported"

    records = []

    for idx, line in enumerate(gb):
        if line.startswith("ACCESSION"):
            access = line[12:].split(" ")[0].strip()


        elif "ORGANISM" in line:
            organism = line[12:].strip()


        elif line.startswith("ORIGIN"):
            seq_idx_start = idx + 1

        elif line.startswith("//"):
            seq_idx_end = idx - 1
            seq_lines_len = seq_idx_end - seq_idx_start
            sequence = []
            for i in range(seq_lines_len+1):
                sequence.append(gb[seq_idx_start + i][10:].strip())
            sequence = "".join(sequence)
            sequence = sequence.replace(" ", "").upper()

            record = GenbankRecord(access, organism, sequence)
            records.append(record)

    return records

def sort_tuples(tuple, p, rev = True):
    """

    :param tuple: List of tuples of 1 or more elements
    :param p: Index of the element that the tuple should be sorted by
    :return: List of tuples sorted by the values of element p of the tuple.
    """
    sorted_tuples = sorted(tuple, key=lambda tup: tup[p], reverse=rev)

    return sorted_tuples

class GenbankRecord:
    """
    Contains information in a GenBank record about a DNA sequence
    """
    def __init__(self, ID, organism, sequence):
        self.ID = ID
        self.organism = organism
        self.sequence = sequence
        self.length = len(self.sequence)
        self.gc = calc_gc(sequence)

    def info(self):
        print "{}\t{}\t{}\t{}".format(self.ID, self.organism, self.gc,
                                      self.length)

def sort_records_gc(records_list):

if __name__ == "__main__":

    """ Parse the Genbank file and store required data """

    records = parse_genbank(argv[1])

    """ Sort and get list of Accession IDs sorted by GC content """
    for record in records:
        gc =
    gc_sorted = sort_tuples(gc, 1)
    sorted_acc = []
    for tup in gc_sorted:
        sorted_acc.append(tup[0])

    """ Write a FASTA document with the sorted sequences """

    fasta_string = ""
    for acc in sorted_acc:
        fasta_string += ">{} {}\n{}\n".format(acc, orgs[acc], seqs[acc])

    fasta_file = open("seqs.fasta", "w")
    fasta_file.write(fasta_string)
    fasta_file.close()

    """ Write a .txt file with the required stats about each sequence """

    tab_file = open("stats.txt", "w")
    tab_string = ""
    for idx, acc in enumerate(sorted_acc):
        tab_string += "".join([acc.ljust(20), str(orgs[acc]).ljust(35),
                               str(gc[idx][1]).ljust(10),
                               str(lengths[acc]).ljust(15), "\n"])

    tab_file.write(tab_string)
    tab_file.close()

    print tab_string


    seq1 = GenbankRecord("BetaClass", "Alien", "CGCGGGTATATAGC")
    seq1.info()

