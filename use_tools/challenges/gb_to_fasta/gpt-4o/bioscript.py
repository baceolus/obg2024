import sys
from Bio import SeqIO

input_file = "input_1.gb"
output_file = "output.txt"

with open(input_file, "r") as gb_file, open(output_file, "w") as fasta_file:
    SeqIO.convert(gb_file, "genbank", fasta_file, "fasta")