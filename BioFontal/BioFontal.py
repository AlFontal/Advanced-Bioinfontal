#!/usr/bin/env python
"""
Script containing some summary Functions
"""

__author__ = "Alejandro Fontal"
__email__ = "jandrofontal@gmail.com"


def calc_gc(sequence):
    """
    :param sequence: String containing a DNA sequence. Case insensitive
    :return: Percentage of GC content rounded with 2 decimals.
    """
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    gc_perc = round(gc_count / len(sequence) * 100, 2)
    return gc_perc


def sort_tuples(tuple, p, rev=True):
    """
    :param tuple: List of tuples of 1 or more elements
    :param p: Index of the element that the tuple should be sorted by
    :return: List of tuples sorted by the values of element p of the tuple.
    """
    sorted_tuples = sorted(tuple, key=lambda tup: tup[p], reverse=rev)
    return sorted_tuples


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


def parse_fastq_file(filename):
    """
    Parses sequences and quality scores from fastQ files

    :param filename: Filename of a fastQ file with Illumina 1.5+ Encoding
    :return: A dictionary with sequences as values and quality scores ranging
    from 0 to 41.
    """
    with open(filename, "r") as file:
        file_lines = filename.readlines()
        seqs_dict = {}
        for idx, line in enumerate(file_lines):
            if line.startswith("@"):
                curr_seq = file_lines[idx + 1].strip()
                seqs_dict[curr_seq] = []
                qc_vals = convert_scores(file_lines[idx + 3].strip())
                seqs_dict[curr_seq] = qc_vals

    return seqs_dict


def run_tool(arg1, arg2, out_fn):
    """
    Runs the tool program in the Linux Shell
    :return: Filename of the output
    """
    path = os.getcwd() + "/" + out_fn
    if os.path.exists(path):
        return out_fn
    else:
        cmd = "tool {} {}".format(arg1, arg2)
        try:
            subprocess.check_call(cmd, shell=True)
            return out_fn
        except subprocess.CalledProcessError:
            print "There was an error - tool command exited with non-zero code"
