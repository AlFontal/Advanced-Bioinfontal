#!/usr/local/bin/python

from __future__ import division
import subprocess
import os
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


def tuple_to_tuple_len(tuples):
    """
    Takes a tuple or a list of tuples and converts it to a tuple or a list of
    tuples that contains the same element as the original in the first
    position and the length of the element in the second position in the
    original.
    :param tuples: tuple or list of tuples
    :return: tuple or list of tuples.
    """
    if type(tuples) == tuple:
        new_tuples = (tuples[0], len(tuples[1]))
        return new_tuples

    elif type(tuples) == list:
        new_tuples = map(lambda tup: (tup[0], len(tup[1])), tuples)
        return new_tuples
    else:
        print "Type Error: Input type not tuples nor list or tuples"


def sort_tuples(tuple, p, rev = True):
    """
    :param tuple: List of tuples of 1 or more elements
    :param p: Index of the element that the tuple should be sorted by
    :return: List of tuples sorted by the values of element p of the tuple.
    """
    sorted_tuples = sorted(tuple, key=lambda tup: tup[p], reverse=rev)

    return sorted_tuples


def run_lastz(ref, query, out_fn="outlastz.txt", format="general"):
    path = os.getcwd() + "/" + out_fn
    if os.path.exists(path):
        return out_fn
    else:
        cmd = "lastz {} {} --output={} --format={}".format(ref, query,
                                                           out_fn, format)
        try:
            subprocess.check_call(cmd, shell=True)
            return out_fn

        except subprocess.CalledProcessError:
           print "There was an error - lastz command exited with non-zero code"


def parse_lastz_out(lastz_fn):
    with open(lastz_fn) as gen_align:
        covered = []
        for line in gen_align.readlines()[1:]:
            covered.append(map(int,line.split("\t")[4:6]))

        covered = sorted(covered, reverse=False)


        """
        for idx, fragment in enumerate(covered):
            if idx == 0:
                continue
            if fragment[1] < covered[idx-1][1]:
                print fragment
        """

class Assembly:
    """
    Object that stores information about assemblies and eases the process of
    storing and retrieving data from them

    """
    def __init__(self, dict):
        """

        :param dict: dictionary containing a parsed assembly file. It should
        contain labels as keys and contigs as values.
        """
        self.sorted_contigs = sort_tuples(tuple_to_tuple_len(dict.items()), 1)
        self.nr_contigs = len(self.sorted_contigs)
        size = 0
        for contig in self.sorted_contigs:
            size += contig[1]
        self.size = size

        current_size = 0
        for idx, contig in enumerate(self.sorted_contigs):
            current_size += contig[1]
            if current_size > self.size/2:
                self.n50_size = contig[1]
                self.n50_idx = idx + 1
                break

    def print_stats(self):
        """
        Prints stats about the assembly
        """
        print "This assembly contains {} contigs and has a total size of {}" \
              " nt, The N50 index is {} and the N50 size {} nt.".format(
            self.nr_contigs, self.size, self.n50_idx, self.n50_size)


if __name__ == "__main__":

    chrom_fn = "chr3.fa"
    contigs_fn = "velvet_15.fa"
    chrom = parse_fasta("chr3.fa")
    contigs = parse_fasta("velvet_15.fa")
    velv = Assembly(contigs)
    velv.print_stats()
    align_output = run_lastz(chrom_fn, contigs_fn)
    parse_lastz_out(align_output)
    a = set(range(1, 4))
    b = set(range(5, 9))
    c = set(range(3, 7))
    d = a | b | c
    print d