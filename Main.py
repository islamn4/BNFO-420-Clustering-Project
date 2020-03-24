"""
from Bio import SeqIO
import timeit
import functions as fun

#non_ambiguous_seq_list = []


for record in SeqIO.parse("fasta_files/ambiguous_aa_test.fasta", "fasta"):
    if(fun.ambiguous_aa(record.seq)):
        SeqIO.write(record.seq, "non_ambiguous_sequences.fasta", "fasta")
    else:
        print(record.seq)
"""

import timeit
from Bio import SeqIO

def ambiguous_aa(seq):
    """
    A function that looks for ambiguous amino acids.
    The ambiguous amino acids include, without the quotes, the following: "B", "Z", "X", "J", "O", "U", ".", "*"

    :parameter
    seq : string
        Contains the amino acid sequence to be searched
    :return:
    True if no ambiguous amino acids are found.
    False if at least one ambiguous amino acid is found.
    """

    seq = seq.upper()

    if(seq.find("B") != -1):
        return False
    elif(seq.find("Z") != -1):
        return False
    elif(seq.find("X") != -1):
        return False
    elif(seq.find("J") != -1):
        return False
    elif(seq.find("O") != -1):
        return False
    elif (seq.find("U") != -1):
        return False
    elif(seq.find(".") != -1):
        return False
    elif(seq.find("*") != -1):
        return False
    else:
        return True

def main_fun():
    for record in SeqIO.parse("fasta_files/ambiguous_aa_test.fasta", "fasta"):
        ambiguous_aa(record.seq)
        #if(ambiguous_aa(record.seq)):
            #SeqIO.write(record.seq, "non_ambiguous_sequences.fasta", "fasta")
        #else:
            #print(record.seq)

seq_file = open("fasta_files/ambiguous_aa_test.fasta", "r")
seq = seq_file.read()
seq_file.close()
seq = seq.strip()

t_aa_fun = timeit.Timer("ambiguous_aa(seq)", setup="from __main__ import seq, ambiguous_aa")
t_main_fun = timeit.Timer("main_fun()", setup="from __main__ import ambiguous_aa, main_fun")

print(t_aa_fun.timeit(100000))
print(t_main_fun.timeit(100000))






























