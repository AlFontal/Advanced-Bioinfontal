#!/usr/bin/env python

"""
Author: Alejandro Fontal
Student Registration Number: 920110-242-090
"""

from P2 import *
from sys import argv
import subprocess
from Bio import Entrez
import os

Entrez.email = "alejandro.fontal@wur.nl"


def wget(link):

    """
    Downloads the content of the provided link to the current folder, basically
    calls the wget function in the unix shell.
    """
    cmd = "wget {}".format(link)

    subprocess.check_output(cmd, shell=True)


def parse_ids(filename):

    names_file = open(filename, "r")
    access = names_file.read().split()

    return access


if __name__ == "__main__":

    path = "{}/p4input.txt".format(os.getcwd())
    if not os.path.exists(path):
        wget("http://www.bioinformatics.nl/courses/BIF-30806/docs/p4input.txt")

    access = parse_ids("p4input.txt")

    handle = Entrez.efetch(db="nucleotide", id=access, rettype="gb",
                               retmode="text")

    genbank_file = handle.readlines()

    stats_list = parse_genbank(genbank_file)

    orgs = stats_list[0]
    access = orgs.keys()
    seqs = stats_list[1]
    lengths = stats_list[3]

    len_tuples = lengths.items()
    shortest_seq = sort_tuples(len_tuples, 1, rev=False)[0][0]

    for key in access:

        print "\n{}\t{}\t{}".format(str(key).ljust(len(max(access)) + 1),
                            str(orgs[key].ljust(len(max(orgs.values())) + 4)),
                            str(lengths[key]).ljust(15))

    """Write FASTA file with the shortest sequence"""

    with open("shortest_seq.fasta", "w") as fasta_file:
        fasta_string = ">{}\t{}\n{}\n".format(shortest_seq,
                                              orgs[shortest_seq],
                                              seqs[shortest_seq])
        print fasta_string
        fasta_file.write(fasta_string)






