import primer3
import sys
from Bio import SeqIO
import json

# Load the test dataset (FASTA file)
fasta_file = "primer_eval/Homo_sapiens_NELFB_sequence.fa"

def design_primers(fasta_file):
    result = {}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(seq_record.seq)
        primers = primer3.designPrimers(
            {
                'SEQUENCE_ID': seq_record.id,
                'SEQUENCE_TEMPLATE': seq
            },
            {
                'PRIMER_MAX_TM': 63.0,
                'PRIMER_MIN_GC': 50.0,
                'PRIMER_NUM_RETURN': 5
            }
        )
        result[seq_record.id] = primers
    return result

# Execute the primer design function, and write to designated output
result = design_primers(fasta_file)
output_path = sys.argv[1]
with open(output_path, "w+") as f:
    json.dump(result, f)