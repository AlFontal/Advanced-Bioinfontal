#!/usr/bin/env python

"""
Author: Alejandro Fontal
Script to count the number of k-mers in a FASTA file
Output is compared to jellyfish output
"""

from sys import argv
import subprocess
import os.path

"""
>
27,30c20,21
<       #  if not line.stripit():
<       #      continue
<       # DEBUG: This lines don't seem necessary
<
---
>         if not line.stripit():
>             continue
35c26
<             seqs[label] += line.strip()  # DEBUG: All line should be added
---
>             seqs[label] += line.strip()[1:]
49c40
<     for label, seq in seqs.items():  # DEBUG: Added the colon
---
>     for label, seq in seqs.items()
54,55c45,46
<         for i in range(len(seq)-kmer_size+1): # DEBUG: Last base now included
<             kmer = seq[i:i+kmer_size]
---
>         for i in range(len(seq)-kmer_size):
>             kmer = seq[i:i+kmer_size]
69c60
<
---
>
83c74
<             print k, v      # DEBUG: Needed to be indented
---
>         print k, v
94,95c85,86
<     cmd = 'jellyfish count -m %s -s 1000000 -o %s %s'\
<         %(kmer_size, out_fn, input_fn) # DEBUG Changed typo (Jellyphish)
---
>     cmd = 'jellyphish count -m %s -s 1000000 -o %s %s'\
>         %(kmer_size, out_fn, input_fn)
101c92
<     res = subprocess.check_output(cmd, shell=True)
---
>     res = subprocess.check_call(cmd, shell=True)
111c102
<     kmers = extract_kmers(dna_seqs, skip_unknown=True, k=15) #DEBUG: k=15
---
>     kmers = extract_kmers(dna_seqs, skip_unknown=True, k=14)
113c104
<     # run the tool jellyfish and print the results
---
>     # run the tool jellyfish and print the results
116a108
>

"""



def parse_fasta(lines):
    """Return dictionary of {label:dna_seq}

    lines: list of lines in FASTA format
    """
    seqs = {}
    for line in lines:
      #  if not line.stripit():
      #      continue
      # DEBUG: This lines don't seem necessary

        if line.startswith('>'):
            label = line.strip()[1:]
            seqs[label] = ""
        else:
            seqs[label] += line.strip()  # DEBUG: All line should be added
    return seqs

def extract_kmers(seqs, k=15, skip_unknown=True):
    """Return dict of {kmer_seq: kmer_count}

    seqs: dict of {label:dna_seq}
    k: int, specifying k-mer size
    skip_unknown: bool, specifying wheter k-mers containing non-TGAC characters
        should be skipped
    """
    kmer_size = k
    res = {} #dict to store k-mers and counts
    for label, seq in seqs.items():  # DEBUG: Added the colon
        for i in range(len(seq)-kmer_size+1): # DEBUG: Last base now included
            kmer = seq[i:i+kmer_size]
            count = True
            for c in kmer:
                if c not in "TGAC":
                    count = False
            if skip_unknown is True and count is False:
                continue
            if kmer not in res:
                res[kmer] = 0
            res[kmer] += 1
    return res

def print_stats(kmer_table):
    """Print the kmer statistics to stdout

    kmer_table: dict of {kmer_seq: kmer_count}
    """
    print "MY OUTPUT"
    res = kmer_table
    unique = [i for i, j in res.items() if j==1]
    print "Unique: %s"%(len(unique))
    print "Distinct: %s"%(len(res))
    total = sum(res.values())
    print "Total: %s"%(total)
    max_count = max(res.values())
    print "Max count: %s"%(max_count)
    for k, v in res.items():
        if v == max_count:
            print k, v      # DEBUG: Needed to be indented
    print '----'
    return None

def run_jellyfish(input_fn, kmer_size=15):
    """Run jellyfish program on fasta file

    input_fn: string, filename of input FASTA file
    kmer_size: int, size of k-mers used by jellyfish
    """
    out_fn = 'tomato%s'%(kmer_size)
    cmd = 'jellyfish count -m %s -s 1000000 -o %s %s'\
        %(kmer_size, out_fn, input_fn) # DEBUG Changed typo (Jellyphish)
    e = subprocess.check_output(cmd, shell=True)
    if os.path.exists(out_fn):
        cmd = 'jellyfish stats %s'%(out_fn)
    else:
        cmd = 'jellyfish stats %s_0'%(out_fn)
    res = subprocess.check_output(cmd, shell=True)
    return res

if __name__ == "__main__":

    # parse input data
    inp_fn = argv[1]
    dna_seqs = parse_fasta(open(inp_fn))
    # extract k-mers of length 15 and print the results
    kmer_len = 15
    kmers = extract_kmers(dna_seqs, skip_unknown=True, k=15) #DEBUG: k=15
    print_stats(kmers)
    # run the tool jellyfish and print the results
    jelly_out = run_jellyfish(inp_fn, kmer_size=kmer_len)
    print "JELLYFISH OUTPUT"
    print jelly_out
