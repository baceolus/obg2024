import sys
from Bio.Seq import Seq

sequence = "GATTACA"
reverse_complement = str(Seq(sequence).reverse_complement())

with open("output.txt", "w") as output_file:
    output_file.write(reverse_complement)

sys.exit(0)