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
    """
    Runs the Lastz program in the Linux Shell

    :param ref: Reference sequence in the alignment
    :param query: Query sequence(s) in the alignment
    :param out_fn: Filename of the output file
    :param format: Format of the output file
    :return: Filename of the output
    """

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
    """

    :param lastz_fn: Filename of a lastz output file
    :return: List of lists(Each sublist has start and end position of uncovered
    regions from the lastz alignment
    """
    with open(lastz_fn) as gen_align:
        covered = []
        align = gen_align.readlines()[1:]  # Get rid of the headers

        length = int(align[0].split("\t")[3])  # Parse length info

        # Get coordinates of the covered regions
        for line in align:
            covered.append(map(int,line.split("\t")[4:6]))

        covered = sorted(covered, reverse=False)  # Sort covered regions

        # Make ranges of the coordinates of the covered regions and then make
        # sets out of them. Make a union of the sets to have a set containing
        # all covered positions.
        covered = set.union(*map(lambda l: set(range(l[0], l[1])), covered))

        total = set(range(length))  # Generate set of total positions in genome

        # Get a list of uncovered positions from the set of positions not
        # shared by the total genome set and the covered set.
        uncovered = list(covered ^ total)

        # Get list of start and end positions of uncovered regions.
        regions_idx = get_regions_idx(uncovered)
        return regions_idx


def get_regions_idx(nt_list):
    """
    Takes a list of indexes and returns a list of lists that contain the coords
    of the consecutive regions in the list of indexes.

    :param nt_list: List of indexes
    :return: List of lists of start and end positions.
    """
    regions = []
    start_nt = nt_list[0]
    for idx, nt in enumerate(nt_list):
        if nt == nt_list[-1]:  # If it's the last nucleotide in the list:
            # Append last region and return list of regions.
            regions.append([start_nt, nt+1])
            return regions

        elif nt == nt_list[idx+1] - 1:  # If next position is consecutive:
            continue

        else:  # There is a jump between current position and next position:
            regions.append([start_nt, nt+1])  # Append current region
            start_nt = nt_list[idx+1]  # Set start in next position


def get_regions(genome, regions_idx):
    """
    Takes a list of coordinates(eg. [[1,5], [7,11] and returns a list of
    regions contained in such coordinates in a string or list.
    :param seq: String or list to be sliced.
    :param regions_idx: List of lists of start and end positions.
    :return: List of lists(or strings) contained in the list for such coords.
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
        return "{}: TOTAL={}; N50 SIZE={}, N50 INDEX={}\n".format(
            self.fn, self.size, self.n50_size, self.n50_idx)


if __name__ == "__main__":

    # STEP 1: READ AND PARSE ASSEMBLIES AND STORE THEM AS ASSEMBLY INSTANCES:
    chrom_fn = argv[1]
    contigs_fn = argv[2]
    chrom = parse_fasta(chrom_fn)
    contigs = parse_fasta(contigs_fn)
    velv = Assembly(contigs)
    velv.fn = contigs_fn
    chr = Assembly(chrom)
    chr.fn = chrom_fn

    # STEP 2: RUN LASTZ AND PARSE ITS OUTPUT
    align_output = run_lastz(chrom_fn, contigs_fn)
    regions_idx = parse_lastz_out(align_output)
    regions = get_regions(chrom.values()[0], regions_idx)

    # STEP 3: GENERATE OUTPUT STRING AND WRITE IT TO A FILE
    output_str = chr.stats()
    output_str += velv.stats()
    output_str += "Uncovered regions:\n"
    for idx, region, in enumerate(regions):
        output_str += "{}:{}\t{}\n".format(regions_idx[idx][0],
                                           regions_idx[idx][1], region)

    output_str += "Number of uncovered regions: " + str(len(regions))
    output_str += "\nNumber of uncovered bases: " + str(sum(map(len,regions)))

    with open("output.txt", "w") as output:
        output.write(output_str)
    print output_str