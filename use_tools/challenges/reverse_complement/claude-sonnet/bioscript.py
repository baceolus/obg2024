import sys

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_comp = ''.join(complement.get(base, base) for base in reversed(seq))
    return rev_comp

sequence = 'GATTACA'
rev_comp_seq = reverse_complement(sequence)

with open('output.txt', 'w') as f:
    f.write(rev_comp_seq)