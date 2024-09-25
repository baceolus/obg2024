import primer3
from Bio import SeqIO
import sys
import json

# Load the test dataset (FASTA file)
fasta_file = "primer_eval/Homo_sapiens_NELFB_sequence.fa"

def design_primers(fasta_file_path):
    """Compute PCR primers for all sequences in a FASTA file."""
    result = {}
    sequences = SeqIO.parse(fasta_file_path, "fasta")
    for seq_record in sequences:
        seq = str(seq_record.seq)
        primers = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': seq_record.id,
                'SEQUENCE_TEMPLATE': seq
            },
            {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 25,
                'PRIMER_OPT_TM': 60.0,
                'PRIMER_MIN_TM': 57.0,
                'PRIMER_MAX_TM': 63.0,
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
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