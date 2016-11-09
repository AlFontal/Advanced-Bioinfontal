#! /usr/bin/env python

from __future__ import division
from sys import argv
import numpy as np
import subprocess

"""
Author: Alejandro Fontal and Roos Goessen
Script to solve P3
"""


def convert_scores(qual_list, offset=64):
    """

    :param qual_list: List of encoded qualities
    :param offset: Offset to take into account (depends on
    the fasq encoding. Default = 64 (Illumina 1.5+))
    :return: List of converted qualities (from 0 to 41)
    """

    qc_vals = []
    for char in qual_list:
        qual = ord(char) - offset
        qc_vals.append(qual)

    return qc_vals


def parse_fastq_file(filename):
    """

    :param filename: Filename of a fastQ file with Illumina 1.5+ Encoding
    :return: A dictionary with sequences as values and quality scores ranging
    from 0 to 41.
    """
    filename = open(filename, "r")
    file = filename.readlines()

    seqs_dict = {}

    for idx, line in enumerate(file):

        if line.startswith("@"):
            curr_seq = file[idx + 1].strip()
            seqs_dict[curr_seq] = []

            qc_vals = convert_scores(file[idx + 3].strip())

            seqs_dict[curr_seq] = qc_vals

    return seqs_dict


def fastq_seq_len(parsed_fq_dict):
    """

    :param parsed_fq_dict: Dictionary obtained from the
    parse_fastq_file function.
    :return: List containing the max length, min length
    and average length of the sequences inside contained
    in the dictionary
    """

    lengths = map(len, parsed_fq_dict.keys())
    max_length = max(lengths)
    min_length = min(lengths)
    avg_length = np.mean(lengths)

    return [max_length, min_length, avg_length]


def get_quality_pos(parsed_fq_dict):
    """

    :param parsed_fq_dict: dictionary with DNA seqs as keys and
    qualities ranging from 0 to 41 as values
    :return: A dictionary with positions as keys and average quality per
    position as values.
    """

    qual_dict = {}
    for read in parsed_fq_dict:
        quality_list = parsed_fq_dict[read]
        for idx, value in enumerate(quality_list):
            if idx in qual_dict.keys():
                qual_dict[idx].append(value)
            else:
                qual_dict[idx] = [value]

    for i in range(len(qual_dict)):
        qual_dict[i] = round(np.mean(qual_dict[i]), 2)

    return qual_dict


def get_diff_qual(normal_qual_dict, trimmed_qual_dict):
    """

    :param normal_dict: Dictionary containing qualities of the original .fq
    :param trimmed_dict: Dictionary containing qualities of the trimmed .fq
    :return: diff: list of integers with the difference between both lists of
    qualities per position.
    """

    diff = []
    for i in range(len(normal_qual_dict)):
        diff.append(round(trimmed_qual_dict[i] - normal_qual_dict[i], 3))

    return diff


def run_trimmer(filename, threshold=30, q=64, out="trimmed.fq"):
    """

    :param filename: Name or path of the fastq file to trim
    :param threshold: Minimum base quality to set the trimming.
    :param q: FASTQ encoding type, default = 64 (Illumina 1.5+)
    :param out: Filename of the output trimmed file. By default
    it is "trimmed.fq"
    :return: Filename of the generated output file.
    """

    cmd = 'fastq_quality_trimmer -t {} -Q {}' \
          ' -i {} -o {}'.format(threshold, q,
                                filename, out)

    subprocess.check_output(cmd, shell=True)

    return out


if __name__ == "__main__":

    original_file = argv[1]

    # Step 1: Get the trimmed fastq file:
    trimmed_fq = run_trimmer(original_file)

    # Step 2: Parse both the original and the trimmed fastq files.
    parsed_orig_fq = parse_fastq_file(original_file)
    parsed_trim_fq = parse_fastq_file(trimmed_fq)

    # Step 3: Calculate length stats for both files.
    max_len_o, min_len_o, avg_len_o = fastq_seq_len(parsed_orig_fq)
    max_len_t, min_len_t, avg_len_t = fastq_seq_len(parsed_trim_fq)

    # Step 4: Calculate average qualities per position and difference between
    # the two files.
    qual_dict_o = get_quality_pos(parsed_orig_fq)
    qual_dict_t = get_quality_pos(parsed_trim_fq)
    diff = get_diff_qual(qual_dict_o, qual_dict_t)

    # Step 5: Print the output to the console:

    print "\nORIGINAL:   max = {}    min = {}    " \
          "avg = {}".format(max_len_o, min_len_o, avg_len_o)

    print "TRIMMED:    max = {}    min = {}    " \
          " avg = {}".format(max_len_t, min_len_t, avg_len_t)

    print "\nThe average read quality per position is:\n "
    print "{}\t{}\t{}\t{}\n".format("Position".center(15),
                                    "Original".center(15),
                                    "Trimmed".center(15),
                                    "Difference".center(15))
    for i in range(max_len_o):
        print "{}\t{:10.2f}\t{:10.2f}\t{:10.2f}".format(str(i + 1).center(15),
                                                qual_dict_o[i],
                                                qual_dict_t[i],
                                                diff[i])

