#!/usr/bin/env python

"""
Author: Alejandro Fontal
Student Registration Number: 920110-242-090
"""

from P2 import *
from sys import argv
import subprocess
from Bio import Entrez

Entrez.email = "alejandro.fontal@wur.nl"


def wget(link):
    """
    Downloads the content of the provided link to the current folder
    """
    cmd = "wget {}".format(link)

    subprocess.check_output(cmd, shell=True)


def parse_ids(filename):

    names_file = open(filename, "r")
    access = names_file.read().split()

    return access

if __name__ == "__main__":


    #wget("http://www.bioinformatics.nl/courses/BIF-30806/docs/p4input.txt")

    access =  parse_ids("p4input.txt")


    parsed_list = []


    for accession_nr in access:
        handle = Entrez.efetch(db="nucleotide", id=accession_nr, rettype="gb",
                               retmode="text")

        gb = handle.readlines()


        parsed_list.append(parse_genbank(gb))

    access = []
    for i in range(len(parsed_list)):
        access.append(parsed_list[i][0].keys()[0])

    lengths = []
    for i in range(len(parsed_list)):
        lengths.append(parsed_list[i][3].items()[0])

    orgs = []
    for i in range(len(parsed_list)):
        orgs.append(parsed_list[i][0].values()[0])

    seqs = []
    for i in range(len(parsed_list)):
        seqs.append(parsed_list[i][1].values()[0])



    shortest_seq = sort_tuples((lengths), 1, rev=False)[0][0]
    idx = access.index(shortest_seq)
    """ Write FASTA file with the shortest sequence """

    with open("shortest_seq.fasta", "w") as fasta_file:
        fasta_string = ">{}\t{}\n{}\n".format(shortest_seq,
                                              orgs[idx],
                                              seqs[idx])

        fasta_file.write(fasta_string)

    tab_string = ""

    for i in range(len(parsed_list)):

        print "{}\t{}\t{}\n".format(str(access[i].ljust(len(max(access))+1)),
                                    str(orgs[i].ljust(len(max(orgs))+4)),
                                    str(lengths[i][1]).ljust(15))








