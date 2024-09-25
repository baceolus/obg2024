import sys

seq = "GATTACA"
rev_comp = "".join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in seq[::-1]])

with open("output.txt", "w") as f:
    f.write(rev_comp)