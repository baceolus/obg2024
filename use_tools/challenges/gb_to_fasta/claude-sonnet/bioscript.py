import sys
from Bio import SeqIO

input_file = "input_1.gb"
output_file = "output.txt"

with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "genbank")
    count = SeqIO.write(sequences, output_handle, "fasta")

print(f"{count} sequences were converted and saved to {output_file}")