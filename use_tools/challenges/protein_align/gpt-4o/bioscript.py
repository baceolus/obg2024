import sys
from Bio import SeqIO, pairwise2

# Read the input FASTA file
input_file = "input_1.fasta"
sequences = list(SeqIO.parse(input_file, "fasta"))

# Extract the first sequence (sequence 1)
seq1 = sequences[0]

# Initialize variables to track the best alignment
best_score = -1
best_id = None

# Iterate over the rest of the sequences
for seq in sequences[1:]:
    # Perform pairwise alignment
    alignments = pairwise2.align.globalxx(seq1.seq, seq.seq)
    # Get the best alignment score
    score = alignments[0][2]
    # Update the best score and ID if this alignment is better
    if score > best_score:
        best_score = score
        best_id = seq.id

# Write the best ID to the output file
with open("output.txt", "w") as output_file:
    output_file.write(best_id + "\n")