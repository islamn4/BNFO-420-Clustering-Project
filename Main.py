"""
3/25
    Currently, this code removes duplicate strains, sequences with ambiguous amino acids, and all strains that are not H5N1.
To-Do:
    1) Counter for each location
    2) Check for duplicate sequences using the strain's accession number
    3) Continue Testing the code! (minimal testing has been conducted so far)
"""
import os
from Bio import SeqIO
import Functions as fun

directory = os.fsencode("not_filtered_fasta_files")
garbage_file = open("test_files/_garbage_sequences.fasta", "w+")
output_file = open("filtered_fasta_files/all_filtered_sequences.fasta", "w+")

human_strain = False
avian_strain = False
swine_strain = False
non_ambiguous_seq_list = []
counter_total_strains = 0
counter_keep = 0
counter_garbage = 0
counter_ambig = 0
counter_duplicate = 0
counter_human = 0
counter_avian = 0
counter_swine = 0

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        output_handle = open("filtered_fasta_files/filtered_" + str(filename) + ".fasta", "w+")

        if(str(filename).upper().find("HUMAN") != -1):
            human_strain = True
        elif(str(filename).upper().find("AVIAN") != -1):
            avian_strain = True
        else:
            swine_strain = True

        for record in SeqIO.parse("not_filtered_fasta_files/" + str(filename), "fasta"):
            counter_total_strains += 1
            if(fun.ambiguous_aa(str(record.seq)) and not(str(record.seq) in non_ambiguous_seq_list) and fun.H5N1_subtype_check(str(record.description))):
                non_ambiguous_seq_list.append(str(record.seq))
                output_handle.write(str(record.__format__("fasta")) + "\n")
                output_file.write(str(record.__format__("fasta")) + "\n")
                if(human_strain):
                    counter_human += 1
                elif(avian_strain):
                    counter_avian += 1
                else:
                    counter_swine += 1
                counter_keep += 1
            else:
                garbage_file.write(str(record.__format__("fasta")) + "\n")
                counter_garbage += 1
                if(not(str(record.seq) in non_ambiguous_seq_list)):
                    counter_duplicate += 1
                else:
                    counter_ambig += 1
        human_strain = False
        avian_strain = False
        swine_strain = False

    else:
        continue
    output_handle.close()
output_file.close()

print("Total number of H5N1 virus strains = " + str(counter_keep))
print("Number of human infecting H5N1 viruses = " + str(counter_human))
print("Number of avian infecting H5N1 viruses = " + str(counter_avian))
print("Number of swine infecting H5N1 viruses = " + str(counter_swine))
print()
print()
print("counter_total_strains = " + str(counter_total_strains))
print("counter_keep + counter_garbage = " + str(counter_keep + counter_garbage))
print("counter_keep = " + str(counter_keep))
print("counter_garbage = " + str(counter_garbage))
print("counter_ambig + counter_duplicates = " + str(counter_ambig + counter_duplicate))
print("counter_ambig = " + str(counter_ambig))
print("counter_duplicates = " + str(counter_duplicate))
















































