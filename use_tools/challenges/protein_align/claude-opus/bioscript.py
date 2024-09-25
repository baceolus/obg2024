#!/usr/bin/env python

from Bio import SeqIO
from Bio.Align import PairwiseAligner

def main():
    # Read input sequences
    sequences = list(SeqIO.parse("input_1.fasta", "fasta"))

    # Initialize best alignment score and sequence id
    best_score = 0
    best_id = None

    # Create pairwise aligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    # Align each sequence to sequence 1
    seq1 = sequences[0]
    for seq2 in sequences[1:]:
        score = aligner.score(seq1.seq, seq2.seq)
        if score > best_score:
            best_score = score
            best_id = seq2.id

    # Write output
    with open("output.txt", "w") as f:
        f.write(best_id + "\n")

if __name__ == "__main__":
    main()