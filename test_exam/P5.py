#!/usr/bin/env python

"""
Author: Alejandro Fontal
Student Registration Number: 920110-242-090
"""

from __future__ import division
from sys import argv
import os.path
import subprocess
import re


def wget(link):

    """
    Downloads the content of the provided link to the current folder, basically
    calls the wget function in the unix shell. Returns the filename of the
    downloaded document as a string.
    """
    cmd = "wget {}".format(link)

    subprocess.check_call(cmd, shell=True)


def fasta_parse(fasta_fn):
    """

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


def calc_ham_dist(seq1, seq2):
    """
    Calculates Hamming distance between two strings of same length.
    :param seq1: String 1
    :param seq2: String 2
    :return: Hamming distance, number of mismatched symbols.
    """
    ham_dist = 0
    for idx, character in enumerate(seq1):
        if character != seq2[idx]:
            ham_dist += 1

    return ham_dist


def calc_perc_identity(seq1, seq2):
    """
    Calculates percentage of equal symbols between two strings of same length.
    :param seq1: String 1
    :param seq2: String 2
    :return: Float with two decimals containing the percentage of identity.
    """
    identities = 0
    for idx, character in enumerate(seq1):
        if character == seq2[idx]:
            identities += 1

    perc_ident = round(((identities/len(seq1))*100),2)
    return perc_ident


def needle(seq_file1, seq_file2,
           gap_open=8, gap_ext=0.5,
           out_file="out.needle"):
    """

    :param seq_file1: File containing DNA sequence to align
    :param seq_file2: File containing DNA sequence(s) to align with the 1st.
    :param gap_open: Gap open penalty to apply during the alignment
    :param gap_ext: Gap extend penalty to apply during the alignment
    :param out_file: Filename of the alignment file generated
    :return: Filename of the out_file
    """
    path = "{}/{}".format(os.getcwd(), out_file)
    if not os.path.exists(path):
        cmd = "needle {} {} -gapopen {} -gapextend {} -outfile {} ".format(
                            seq_file1, seq_file2, gap_open, gap_ext, out_file)

        status = subprocess.check_call(cmd, shell=True)
        print "Exit status {} for Needle".format(status)
    return out_file

def split_records(align_fn):
    """

    :param align_fn: Filename of a a file containing the output of Needle.
    :return: List of records (string)
    """
    records = []
    record = ""
    with open(align_fn) as align_file:
        align_list = align_file.readlines()
        idx = 0
        while idx < len(align_list)-1:
            if align_list[idx].startswith("# Aligned_sequences:"):
                if len(records) > 0:
                    records.append(record)
                    record = ""
                    idx += 1
                else:
                    records.append(record)
                    record = "".join(align_list[idx-2:idx+1])
                    idx += 1
            else:
                record += align_list[idx]
                if align_list[idx].startswith("#----------------"
                                              "-----------------------"):
                    records.append(record)
                    return records[1:]
                idx += 1


def extract_alignment(align_record):
    """

    :param align_record: String containing a record of the output of an
    alignment made with needle.
    :return: String containing the symbols between the sequences aligned.
    """
    align_str = ""
    matches = re.findall(r"[\s]{21}[|: .]{1,50}[\n]", align_record)

    for match in matches:
        align_str += match[21:len(match)-1]

    return align_str


def get_labels_record(align_record):
    """

    :param align_record: String containing one record of the output of an
    alignment made with needle.
    :return: List of labels(strings) of the sequences in the record.
    """

    matches = re.findall(r"#\s\d:\s[\w]{4,15}", align_record)

    labels = []
    for text in matches:
        label = text[5:]
        labels.append(label)

    return labels

def get_aligned_seqs(record):
    """

    :param record: String containing one record of the output of an
    alignment made with needle.
    :return: List containing the two aligned sequences.
    """

    seqs = []
    for seq in get_labels_record(record):
        regexp = re.escape(seq) + r'[\s]{5,15}[\d]{1,3}[\s][\w|-]{1,50}'
        matches = re.findall(regexp, record)
        seq = ""
        for match in matches:
            seq += match[21:]
        seqs.append(seq)

    return seqs


def get_seqs_from_align_record(parsed_dict, record):
    """

    :param parsed_dict: Dictionary containing keys as labels and sequences as
    values
    :param record: String containing one record of the output of an
    alignment made with needle.
    :return: list with the original sequences aligned in the record.
    """

    keys = []
    seqs = []
    for seq in get_labels_record(record):
        for key in parsed_total.keys():
            if seq in key:
                keys.append(key)
                break

    for key in keys:
        seqs.append(parsed_dict[key])

    return seqs



if __name__ == "__main__":

    """
    In case files haven't been downloaded, check if required files are on
    current working directory, if they are not, download them.
    """

    path = "{}/ref.fasta".format(os.getcwd())
    if not os.path.exists(path):
        wget("http://www.bioinformatics.nl/"
             "courses/BIF-30806/docs/ref.fasta")

    path2 = "{}/related.fasta".format(os.getcwd())
    if not os.path.exists(path2):
        wget("http://www.bioinformatics.nl/"
             "courses/BIF-30806/docs/related.fasta")


    ref_fn = argv[1]
    rel_fn = argv[2]

    # Parse both fasta files and store seqs in dictionaries.
    parsed_ref = fasta_parse(ref_fn)
    parsed_rel = fasta_parse(rel_fn)

    # Merge both dictionaries in one
    parsed_ref.update(parsed_rel)
    parsed_total = parsed_ref

    # Run needle from shell to generate the alignment file
    align_file = needle(ref_fn, rel_fn)

    # Split the alignment file in a record per record basis and store in list.
    records = split_records(align_file)


    print "\n{}\t{}\t{}\t{}\t{}\t{}\n".format("Sequence1".center(15),
                                            "Length".center(15),
                                            "Sequence2".center(15),
                                            "Length".center(15),
                                            "Hamming".center(15),
                                            "Identity".center(15))
    print "-------------------------------"*3 +"\n"

    for record in records:  # For every alignment in the alignment file:

        label_list = get_labels_record(record)  # Get labels of seqs aligned

        # Get the original sequences of the seqs aligned in the record.
        seqs = get_seqs_from_align_record(parsed_total, record)

        # Get the aligned seqs to perform calculations on them
        aligned_seqs = get_aligned_seqs(record)

        # Calculate % of identity between seqs
        identity = calc_perc_identity(aligned_seqs[0], aligned_seqs[1])

        # Calculate Hamming distance between seqs
        hamm = calc_ham_dist(aligned_seqs[0], aligned_seqs[1])

        print "{}\t{}\t{}\t{}\t{}\t{}\n".format(label_list[0].center(15),
                                            str(len(seqs[0])).center(15),
                                            label_list[1].center(15),
                                            str(len(seqs[1])).center(15),
                                            str(hamm).center(15),
                                            str(identity).center(15))












