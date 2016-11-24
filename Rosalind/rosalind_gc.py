"""
Author: Alejandro Fontal
"""


from __future__ import division
from sys import argv
fileName = argv[1]


def fasta_parse(fileName):
    """

    :param fileName: Name of a FASTA file
    :return: Two lists, one containing names of the seqs and one with the seqs.
    """

    fastaFile = open(fileName, 'r')
    seq = []
    seqs = []
    names = []
    for line in fastaFile:
        if line[0] == ">":
            line = line.rstrip()
            names.append(line[1:])
            if seq:
                seq = "".join(seq)
                seqs.append(seq)
                seq = []
        else:
            line = line.rstrip()
            seq.append(line)
    seq = "".join(seq)
    seqs.append(seq)
    fastaFile.close()
    return names, seqs

def calc_max_GC(seqs, names):
    """
    Takes a list of sequences and gives the GC% of the one with the most GC%

    :param seqs: List of DNA sequences in string form
    :param names: List of DNA sequence names.
    :return: Two values. The maximum percentage of GC content of the sequences
    and the name of the sequence with the maximum percentage.
    """

    gc_contents = {}
    for idx, seq in enumerate(seqs):
        gc = seq.count("G") + seq.count("C")
        gc_perc = (gc/len(seq)) * 100
        gc_contents[gc_perc] = names[idx]

    max_gc = max(gc_contents.keys())
    max_seq = gc_contents[max_gc]

    return max_gc, max_seq

names, seqs = fasta_parse(fileName)

max_gc, max_seq = calc_max_GC(seqs, names)
f = open('output.txt', 'w')
f.write(str(max_seq))
f.write("\n")
f.write(str(max_gc))
f.close()
print str(max_seq)
print max_gc