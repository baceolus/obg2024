import sys
from Bio import SeqIO

# Read the GenBank file and convert to FASTA format
with open("input_1.gb", "r") as input_handle, open("output.txt", "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "genbank")
    SeqIO.write(sequences, output_handle, "fasta")

print("Conversion complete. FASTA file saved as output.txt")