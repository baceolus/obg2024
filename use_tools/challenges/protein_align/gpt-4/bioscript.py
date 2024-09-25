```python
import os
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

# Read the fasta file
sequences = list(SeqIO.parse("input_1.fasta", "fasta"))

# Get the first sequence
seq1 = sequences[0].seq

# Initialize the best alignment score and sequence ID
best_score = 0
best_id = ""

# Iterate over the rest of the sequences
for seq_record in sequences[1:]:
    # Perform pairwise alignment
    alignments = pairwise2.align.globalds(seq1, seq_record.seq, matlist.blosum62, -10, -0.5)
    
    # Get the score of the best alignment
    score = alignments[0].score
    
    # If this score is better than the current best, update the best score and sequence ID
    if score > best_score:
        best_score = score
        best_id = seq_record.id

# Write the ID of the sequence with the best alignment to the output file
with open("output.txt", "w") as output_file:
    output_file.write(best_id)
```