import Bio.Seq

sequence = Bio.Seq.Seq("GATTACA")
reverse_complement = sequence.reverse_complement()

with open("output.txt", "w") as output_file:
    output_file.write(str(reverse_complement))