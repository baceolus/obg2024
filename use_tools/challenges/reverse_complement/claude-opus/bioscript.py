import os

seq = "GATTACA"
rev_comp = seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

with open("output.txt", "w") as outfile:
    outfile.write(rev_comp)