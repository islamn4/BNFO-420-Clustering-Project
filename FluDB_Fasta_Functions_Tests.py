# Author: Stephen Shea
# Created: 4/25/20

import os
from Strain import Strain
import NCBI_Fasta_Functions as ncbi
import FluDB_Fasta_Functions as flu

# folder = os.fsencode("Fasta_Files")
folder = os.fsencode("Fasta_Files_Test")
seasonal_folder = os.fsencode("Fasta_Files/seasonal_analysis")

filename_list = os.listdir(folder)
filename_list.sort()
seasonal_filename_list = os.listdir(seasonal_folder)
seasonal_filename_list.sort()

for file in filename_list:

    filename = os.fsdecode(file)
    counter = 1

    if filename.endswith(".fasta") or filename.endswith(".fa"):
        str_file = ""
        with open("Fasta_Files/" + str(filename), "r") as str_file_temp:
            str_file = str_file_temp.read()
            str_file_temp.close()

        for fasta_sequence in str_file.split(">"):
            header = fasta_sequence[:fasta_sequence.find("\n")].strip()
            seq = fasta_sequence[len(header):].strip().upper().replace("\n", "")
            # print("Header:")
            # print(header)
            # print()
            # print("Sequence:")
            # print(seq)
            # print()
            if len(fasta_sequence) < 10:
                # print("*****:" + str(fasta_sequence) + ":****")
                continue
            elif flu.parse_year(header) == None:
                print("Strain " + str(counter) + ":\t\t\t\t\t" + str(flu.parse_year(header)))
                print(header)
            else:
                print("Strain " + str(counter) + ":\t\t\t\t\t" + str(flu.parse_year(header)))
            # print(flu.parse_year(header))
            print()
            counter += 1













