from Bio import SeqIO

def convert_gb_to_fasta(input_file, output_file):
    records = SeqIO.parse(input_file, "genbank")
    SeqIO.write(records, output_file, "fasta")

convert_gb_to_fasta("input_1.gb", "output.txt")