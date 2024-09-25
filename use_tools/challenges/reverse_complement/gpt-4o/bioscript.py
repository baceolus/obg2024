from Bio.Seq import Seq

sequence = "GATTACA"
seq_obj = Seq(sequence)
reverse_complement = str(seq_obj.reverse_complement())

with open("output.txt", "w") as output_file:
    output_file.write(reverse_complement)