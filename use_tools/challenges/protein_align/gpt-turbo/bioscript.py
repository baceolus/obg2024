```python
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

def main():
    input_file = "input_1.fasta"
    output_file = "output.txt"
    
    # Read all sequences from the fasta file
    sequences = list(SeqIO.parse(input_file, "fasta"))
    
    if len(sequences) < 2:
        print("Not enough sequences to perform alignment.")
        sys.exit(1)
    
    # Sequence 1 to align with others
    seq1 = sequences[0].seq
    
    # Initialize aligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    best_score = -1
    best_id = None
    
    # Compare sequence 1 with all other sequences
    for seq_record in sequences[1:]:
        alignment_score = aligner.score(seq1, seq_record.seq)
        if alignment_score > best_score:
            best_score = alignment_score
            best_id = seq_record.id
    
    # Write the best matching sequence ID to the output file
    with open(output_file, "w") as output:
        output.write(best_id)

if __name__ == "__main__":
    main()
```