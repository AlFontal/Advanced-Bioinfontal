# C:/Python27
"""
Author: Alejandro Fontal
Student Registration Number: 920110242090
Script
"""
from sys import argv

fileName = argv[1]
dna_file = open(fileName)
dna_seq = dna_file.read()

def count_nucleotides(dna_seq):
    """
    Counts occurences of each nucleotide in a DNA seq and outputs the count.

    """

    A = dna_seq.count("A")
    C = dna_seq.count("C")
    G = dna_seq.count("G")
    T = dna_seq.count("T")
    return A, C, G, T


A, C, G, T = count_nucleotides(dna_seq)

output_string = '{} {} {} {}'.format(A, C, G, T)
output_dna = open('output_dna.txt', "w")
output_dna.write(output_string)
print 'The output should look like this:\n'
print output_string