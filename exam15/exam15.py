#!/usr/local/bin/python

from __future__ import division
import subprocess
import os
from sys import argv

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
        align = gen_align.readlines()[1:]
        length = int(align[0].split("\t")[3])
        for line in align:
            covered.append(map(int,line.split("\t")[4:6]))

        covered = sorted(covered, reverse=False)


        covered = set.union(*map(lambda l: set(range(l[0], l[1])), covered))

        total = set(range(length))
        uncovered = list(covered ^ total)
        regions_idx = get_regions_idx(uncovered)
        return regions_idx

def get_regions_idx(nt_list):
    regions = []
    start_nt = nt_list[0]
    for idx, nt in enumerate(nt_list):
        if nt == nt_list[-1]:
            regions.append([start_nt, nt+1])
            return regions

        elif nt == nt_list[idx+1] - 1:
            continue

        else:
            regions.append([start_nt, nt+1])
            start_nt = nt_list[idx+1]


def get_regions(genome, regions_idx):
    """

    :param seq:
    :param regions_idx:
    :return:
    """
    seqs = []
    for region in regions_idx:
        seq = genome[region[0]:region[1]]
        seqs.append(seq)

    return seqs


class Assembly:
    """
    Object that stores information about assemblies and eases the process of
    storing and retrieving data from them

    """
    def __init__(self, contigs):
        """

        :param dict: dictionary containing a parsed assembly file. It should
        contain labels as keys and contigs as values.
        """
        self.nr_contigs = len(contigs)
        self.sorted_contigs = \
            sort_tuples(tuple_to_tuple_len(contigs.items()), 1)

        self.size = self.calc_size()
        self.n50_size, self.n50_idx = self.nk(50)

        self.fn = ""

    def calc_size(self):
        size = 0
        for contig in self.sorted_contigs:
            size += contig[1]
        return size

    def nk(self, k):
        """
        """
        current_size = 0
        for idx, contig in enumerate(self.sorted_contigs):
            current_size += contig[1]
            if current_size > self.size * (k / 100):
                nk_size = contig[1]
                nk_index = idx + 1
                return nk_size, nk_index

    def stats(self):
        """
        Prints stats about the assembly
        """
        return "{}: TOTAL={}; N50 SIZE={}, N50 INDEX={}".format(
            self.fn, self.size, self.n50_size, self.n50_idx)


if __name__ == "__main__":

    chrom_fn = argv[1]
    contigs_fn = argv[2]
    chrom = parse_fasta(chrom_fn)
    contigs = parse_fasta(contigs_fn)
    velv = Assembly(contigs)
    velv.fn = contigs_fn
    chr = Assembly(chrom)
    chr.fn = chrom_fn
    align_output = run_lastz(chrom_fn, contigs_fn)
    regions_idx = parse_lastz_out(align_output)
    regions = get_regions(chrom.values()[0], regions_idx)

    print chr.stats()
    print velv.stats()
    print "Uncovered regions:"
    for idx, region, in enumerate(regions):
        print "{}:{}\t{}".format(regions_idx[idx][0], regions_idx[idx][1],
                                 region)

    print "Number of uncovered regions: ", str(len(regions))
    print "Number of uncovered bases: ", str(sum(map(len,regions)))
