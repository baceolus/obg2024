import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Read input sequences
sequences = list(SeqIO.parse("input_1.fasta", "fasta"))

# Get sequence 1 as query
query_seq = sequences[0].seq

# Perform local blastp against all other sequences
result_handle = NCBIWWW.qblast("blastp", "nr", query_seq)

# Parse results
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)

# Find best alignment
best_align = None
best_score = 0
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.score > best_score:
            best_score = hsp.score
            best_align = alignment

# Write output
if best_align:
    with open("output.txt", "w") as f:
        f.write(best_align.title.split("|")[3])
else:
    print("No significant alignments found")