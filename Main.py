"""
3/25
    Currently, this code removes duplicate strains, sequences with ambiguous amino acids, and all strains that are not H5N1.
To-Do:
    *) The
    1) Counter for each location
    2) Check for duplicate sequences using the strain's accession number
    3) Continue Testing the code! (minimal testing has been conducted so far)
    4) Create a different file for each subtype
    5) Make sure duplicates that are on two different input files will still be caught and removed.
    6) Place strains with questionable subtypes(ex: H9) into a dedicated file so the strains can be reviewed by hand.
    7)
"""
""" Code from clustering idea
        import os
        from Bio import SeqIO
        import Functions as fun
        
        
        directory = os.fsencode("not_filtered_fasta_files")
        test_file = open("test_files/test.fasta", "w+")
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
                output_handle = open("filtered_fasta_files/filtered_" + str(filename), "w+")
        
                if(str(filename).upper().find("HUMAN") != -1):
                    human_strain = True
                elif(str(filename).upper().find("AVIAN") != -1):
                    avian_strain = True
                else:
                    swine_strain = True
        
                for record in SeqIO.parse("not_filtered_fasta_files/" + str(filename), "fasta"):
                    counter_total_strains += 1
                    if(fun.ambiguous_aa(str(record.seq)) and not(str(record.seq) in non_ambiguous_seq_list) and (fun.subtype_filter(str(record.description)) != "none")):
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
"""

import os
from Strain import Strain as strain
import NCBI_Fasta_Functions as ncbi
import Functions as fun

folder = os.fsencode("not_filtered_fasta_files")
test_file = open("test_files/test.fasta", "w+")
output_file = open("filtered_fasta_files/all_filtered_sequences.fasta", "w+")

strain_list = []

counter_total_strains = 0

for file in os.listdir(folder):
    human_strain = False
    filename = os.fsdecode(file)
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        output_handle = open("filtered_fasta_files/filtered_" + str(filename), "w+")
        str_file = ""

        if (str(filename).upper().find("HUMAN") != -1):
            human_strain = True
        with open("not_filtered_fasta_files/" + str(filename), "r") as str_file_temp:
            str_file = str_file_temp.read()
            str_file_temp.close()
        for fasta_sequence in str_file.split(">"):
            header = fasta_sequence[:fasta_sequence.find("\n")].strip()
            seq = fasta_sequence[len(header):].strip().upper().replace("\n", "")
            subtype = ncbi.parse_subtype(header)

            if(ncbi.check_ambiguous_bases(seq) or subtype == None):
                continue
            counter_total_strains += 1
            accession_num = ncbi.parse_id_num(header)
            if human_strain:
                host = "Human"
            else: host = "Avian"
            gene = str(filename)

            if(gene[0] == "M"):
                gene = gene[:gene.find("2") + 1]
            elif(gene.find("N") == 0 and gene[1] == "S"):
                gene = gene[:gene.find("2") + 1]
            else: gene = gene[:gene.find("_")]

            country = ncbi.parse_country(header)
            year = ncbi.parse_year(header)
            # print(counter_total_strains)
            # print(header)
            # print(fasta_sequence)
            # print()
            # print(seq)
            # print()
            # print()

            temp_strain = strain(gene, subtype, host, country, year, accession_num, seq)
            strain_list.append(temp_strain)







    output_handle.close()
output_file.close()
#print(strain.CODON_DICT.values())

for i in strain_list: print(i.get_codon_matrix())

























