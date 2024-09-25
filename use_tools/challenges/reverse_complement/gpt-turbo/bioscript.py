import Bio.Seq

sequence = "GATTACA"
rev_comp_sequence = str(Bio.Seq.Seq(sequence).reverse_complement())

with open("output.txt", "w") as file:
    file.write(rev_comp_sequence)