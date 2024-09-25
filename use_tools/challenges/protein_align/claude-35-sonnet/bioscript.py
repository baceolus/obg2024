import sys
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

def align_sequences(seq1, seq2):
    alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
    return alignment.score

def main():
    sequences = list(SeqIO.parse("input_1.fasta", "fasta"))
    
    if len(sequences) < 2:
        print("Error: Not enough sequences in the input file.")
        sys.exit(1)

    reference_seq = sequences[0].seq
    best_score = -1
    best_id = ""

    for seq_record in sequences[1:]:
        score = align_sequences(reference_seq, seq_record.seq)
        if score > best_score:
            best_score = score
            best_id = seq_record.id

    with open("output.txt", "w") as output_file:
        output_file.write(best_id)

if __name__ == "__main__":
    main()