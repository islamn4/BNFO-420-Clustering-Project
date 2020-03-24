
from Bio import SeqIO
import Functions as fun

#non_ambiguous_seq_list = []


for record in SeqIO.parse("ambiguous_aa_test.fasta", "fasta"):
    if(fun.ambiguous_aa(record.seq)):
        SeqIO.write(record.seq, "non_ambiguous_sequences.fasta", "fasta")
    else:
        print(record.seq)

"""""
from Bio import SeqIO
import timeit

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

"""




























