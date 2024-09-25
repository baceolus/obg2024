#!/usr/bin/env python

from Bio import SeqIO

with open("input_1.gb", "r") as input_handle, open("output.txt", "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "genbank")
    SeqIO.write(sequences, output_handle, "fasta")