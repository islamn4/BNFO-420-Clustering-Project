#Author: Stephen Shea
#Created: 3/2020
"""

"""

import os
from Strain import Strain
import NCBI_Fasta_Functions as ncbi
import FluDB_Fasta_Functions as flu
from xlwt import Workbook

folder = os.fsencode("Fasta_Files")
seasonal_folder = os.fsencode("Fasta_Files/seasonal_analysis")


counter_total_strains = 0
counter_no_subtype = 0
counter_ambiguous_strains = 0
counter_duplicate_strains = 0
counter_strains_post_filtering = 0

# strain_list = []

H5N1_list = []
H7N9_list = []
H1N1_list = []
H3N2_list = []
H1N1_07_08_list = []
H1N1_08_09_list = []
H3N2_10_11_list = []
H3N2_11_12_list = []

H5N1_filename_list = []
H7N9_filename_list = []
H1N1_filename_list = []
H3N2_filename_list = []
H1N1_07_08_filename_list = []
H1N1_08_09_filename_list = []
H3N2_10_11_filename_list = []
H3N2_11_12_filename_list = []

filename_list = os.listdir(folder)
filename_list.sort()
seasonal_filename_list = os.listdir(seasonal_folder)
seasonal_filename_list.sort()

for file in filename_list:
    human_strain = False
    avian_strain = False
    swine_strain = False
    H5N1_strain = False
    H7N9_strain = False
    H1N1_strain = False
    H3N2_strain = False
    filename = os.fsdecode(file)

    if filename.endswith(".fasta") or filename.endswith(".fa"):
        print(str(filename))
        str_file = ""
        if str(filename).upper().find("HUMAN") != -1:
            human_strain = True
        elif str(filename).upper().find("AVIAN") != -1:
            avian_strain = True
        elif str(filename).upper().find("SWINE") != -1:
            swine_strain = True

        if str(filename).upper().find("H5N1") != -1:
            H5N1_strain = True
            H5N1_list.append({str(filename):[]})
            H5N1_filename_list.append(str(filename))
        elif str(filename).upper().find("H7N9") != -1:
            H7N9_strain = True
            H7N9_list.append({str(filename):[]})
            H7N9_filename_list.append(str(filename))
        elif str(filename).upper().find("H1N1") != -1:
            H1N1_strain = True
            H1N1_list.append({str(filename):[]})
            H1N1_filename_list.append(str(filename))
        elif str(filename).upper().find("H3N2") != -1:
            H3N2_strain = True
            H3N2_list.append({str(filename):[]})
            H3N2_filename_list.append(str(filename))

        with open("Fasta_Files/" + str(filename), "r") as str_file_temp:
            str_file = str_file_temp.read()
            str_file_temp.close()

        for fasta_sequence in str_file.split(">"):
            if len(fasta_sequence) < 5:
                # print("TEST!!!!!!!!!!!!!!")
                # print(str(fasta_sequence))
                # print("TEST!!!!!!!!!!!!!!")
                continue
            duplicate = False
            counter_total_strains += 1
            header = fasta_sequence[:fasta_sequence.find("\n")].strip()
            seq = fasta_sequence[len(header):].strip().upper().replace("\n", "")

            # Set strain subtype
            subtype = ncbi.parse_subtype(header)

            # Check for ambiguous bases
            if ncbi.check_ambiguous_bases(seq):
                counter_ambiguous_strains += 1
                continue

            # Check for subtype
            if subtype == None:
                counter_no_subtype += 1
                continue
            # Set strain accession number
            accession_num = flu.parse_id_num(header)

            # Set strain host
            if human_strain:
                host = "Human"
            elif avian_strain:
                host = "Avian"
            elif swine_strain:
                host = "Swine"
            else: host = None

            # Set strain gene
            gene = str(filename)
            gene = gene[:gene.find(".")]
            gene = gene[gene.find("_",8)+1:]

            # flu.parse_gene(header)

            # Set strain country
            country = ncbi.parse_country(header)
            # Set strain year
            year = ncbi.parse_year(header)

            temp_strain = Strain(gene, subtype, host, country, year, accession_num, seq)

            """
            if H5N1_strain:
                H5N1_list[len(H5N1_list)-1].get(str(filename)).append(temp_strain)
            elif H7N9_strain:
                H7N9_list[len(H7N9_list)-1].get(str(filename)).append(temp_strain)
            elif H1N1_strain:
                H1N1_list[len(H1N1_list)-1].get(str(filename)).append(temp_strain)
            elif H3N2_strain:
                H3N2_list[len(H3N2_list)-1].get(str(filename)).append(temp_strain)
            else: print("temp_strain object was not added a subtype list")

            strain_list.append(temp_strain)
            """

            # Check for duplicate strains
            if H5N1_strain:
                for strain in H5N1_list[len(H5N1_list) - 1].get(str(filename)):
                    if temp_strain.get_accession_num() == strain.get_accession_num():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                    elif temp_strain.get_seq() == strain.get_seq():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                if not duplicate:
                    H5N1_list[len(H5N1_list) - 1].get(str(filename)).append(temp_strain)
                    # strain_list.append(temp_strain)
            elif H7N9_strain:
                for strain in H7N9_list[len(H7N9_list) - 1].get(str(filename)):
                    if temp_strain.get_accession_num() == strain.get_accession_num():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                    elif temp_strain.get_seq() == strain.get_seq():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                if not duplicate:
                    H7N9_list[len(H7N9_list) - 1].get(str(filename)).append(temp_strain)
                    # strain_list.append(temp_strain)
            elif H1N1_strain:
                for strain in H1N1_list[len(H1N1_list) - 1].get(str(filename)):
                    if temp_strain.get_accession_num() == strain.get_accession_num():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                    elif temp_strain.get_seq() == strain.get_seq():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                if not duplicate:
                    H1N1_list[len(H1N1_list) - 1].get(str(filename)).append(temp_strain)
                    # strain_list.append(temp_strain)
            elif H3N2_strain:
                for strain in H3N2_list[len(H3N2_list) - 1].get(str(filename)):
                    if temp_strain.get_accession_num() == strain.get_accession_num():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                    elif temp_strain.get_seq() == strain.get_seq():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                if not duplicate:
                    H3N2_list[len(H3N2_list) - 1].get(str(filename)).append(temp_strain)
                    # strain_list.append(temp_strain)
            else:
                print("temp_strain object was not added to a subtype list")
                # strain_list.append(temp_strain)
            if not duplicate:
                counter_strains_post_filtering += 1

print()
H5N1_filename_list.sort()
H7N9_filename_list.sort()
H1N1_filename_list.sort()
H3N2_filename_list.sort()

# Organize the 4 subtype lists
count = 1
for i in range(int(len(H5N1_list)/2), len(H5N1_list), 1):
    H5N1_list.insert(count, H5N1_list.pop(i))
    H5N1_filename_list.insert(count, H5N1_filename_list.pop(i))
    count += 2
count = 1
for i in range(int(len(H7N9_list)/2), len(H7N9_list), 1):
    H7N9_list.insert(count, H7N9_list.pop(i))
    H7N9_filename_list.insert(count, H7N9_filename_list.pop(i))
    count += 2
count = 0
for i in range(int(len(H1N1_list)/2), len(H1N1_list), 1):
    H1N1_list.insert(count, H1N1_list.pop(i))
    H1N1_filename_list.insert(count, H1N1_filename_list.pop(i))
    count += 2
count = 0
for i in range(int(len(H3N2_list)/2), len(H3N2_list), 1):
    H3N2_list.insert(count, H3N2_list.pop(i))
    H3N2_filename_list.insert(count, H3N2_filename_list.pop(i))
    count += 2

#
#########################################################################
#
#########################################################################


for file in seasonal_filename_list:
    human_strain = False
    avian_strain = False
    swine_strain = False
    H1N1_07_08_strain = False
    H1N1_08_09_strain = False
    H3N2_10_11_strain = False
    H3N2_11_12_strain = False
    filename = os.fsdecode(file)

    if filename.endswith(".fasta") or filename.endswith(".fa"):
        print(str(filename))
        str_file = ""
        if str(filename).upper().find("HUMAN") != -1:
            human_strain = True
        elif str(filename).upper().find("AVIAN") != -1:
            avian_strain = True
        elif str(filename).upper().find("SWINE") != -1:
            swine_strain = True

        if str(filename).upper().find("07_08") != -1:
            H1N1_07_08_strain = True
            H1N1_07_08_list.append({str(filename): []})
            H1N1_07_08_filename_list.append(str(filename))
        elif str(filename).upper().find("08_09") != -1:
            H1N1_08_09_strain = True
            H1N1_08_09_list.append({str(filename): []})
            H1N1_08_09_filename_list.append(str(filename))
        elif str(filename).upper().find("10_11") != -1:
            H3N2_10_11_strain = True
            H3N2_10_11_list.append({str(filename): []})
            H3N2_10_11_filename_list.append(str(filename))
        elif str(filename).upper().find("11_12") != -1:
            H3N2_11_12_strain = True
            H3N2_11_12_list.append({str(filename): []})
            H3N2_11_12_filename_list.append(str(filename))

        with open("Fasta_Files/seasonal_analysis/" + str(filename), "r") as str_file_temp:
            str_file = str_file_temp.read()
            str_file_temp.close()

        for fasta_sequence in str_file.split(">"):
            if len(fasta_sequence) < 5:
                # print("TEST!!!!!!!!!!!!!!")
                # print(str(fasta_sequence))
                # print("TEST!!!!!!!!!!!!!!")
                continue
            duplicate = False
            counter_total_strains += 1
            header = fasta_sequence[:fasta_sequence.find("\n")].strip()
            seq = fasta_sequence[len(header):].strip().upper().replace("\n", "")

            # Set strain subtype
            subtype = ncbi.parse_subtype(header)

            # Check for ambiguous bases
            if ncbi.check_ambiguous_bases(seq):
                counter_ambiguous_strains += 1
                continue

            # Check for subtype
            if subtype == None:
                counter_no_subtype += 1
                continue
            # Set strain accession number
            accession_num = flu.parse_id_num(header)

            # Set strain host
            if human_strain:
                host = "Human"
            elif avian_strain:
                host = "Avian"
            elif swine_strain:
                host = "Swine"
            else:
                host = None

            # Set strain gene
            gene = str(filename)
            gene = gene[:gene.find(".")]
            gene = gene[gene.find("_", 8) + 1:]
            gene = gene[:gene.find("_")]

            # flu.parse_gene(header)

            # Set strain country
            country = ncbi.parse_country(header)
            # Set strain year
            year = ncbi.parse_year(header)

            temp_strain = Strain(gene, subtype, host, country, year, accession_num, seq)

            # Check for duplicate strains
            if H1N1_07_08_strain:
                for strain in H1N1_07_08_list[len(H1N1_07_08_list) - 1].get(str(filename)):
                    if temp_strain.get_accession_num() == strain.get_accession_num():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                    elif temp_strain.get_seq() == strain.get_seq():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                if not duplicate:
                    H1N1_07_08_list[len(H1N1_07_08_list) - 1].get(str(filename)).append(temp_strain)
                    # strain_list.append(temp_strain)
            elif H1N1_08_09_strain:
                for strain in H1N1_08_09_list[len(H1N1_08_09_list) - 1].get(str(filename)):
                    if temp_strain.get_accession_num() == strain.get_accession_num():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                    elif temp_strain.get_seq() == strain.get_seq():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                if not duplicate:
                    H1N1_08_09_list[len(H1N1_08_09_list) - 1].get(str(filename)).append(temp_strain)
                    # strain_list.append(temp_strain)
            elif H3N2_10_11_strain:
                for strain in H3N2_10_11_list[len(H3N2_10_11_list) - 1].get(str(filename)):
                    if temp_strain.get_accession_num() == strain.get_accession_num():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                    elif temp_strain.get_seq() == strain.get_seq():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                if not duplicate:
                    H3N2_10_11_list[len(H3N2_10_11_list) - 1].get(str(filename)).append(temp_strain)
                    # strain_list.append(temp_strain)
            elif H3N2_11_12_strain:
                for strain in H3N2_11_12_list[len(H3N2_11_12_list) - 1].get(str(filename)):
                    if temp_strain.get_accession_num() == strain.get_accession_num():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                    elif temp_strain.get_seq() == strain.get_seq():
                        counter_duplicate_strains += 1
                        duplicate = True
                        break
                if not duplicate:
                    H3N2_11_12_list[len(H3N2_11_12_list) - 1].get(str(filename)).append(temp_strain)
                    # strain_list.append(temp_strain)
            else:
                print("temp_strain object was not added to a seasonal subtype list")
                # strain_list.append(temp_strain)
            if not duplicate:
                counter_strains_post_filtering += 1

print()
#
# for i in range(0, int(len(H1N1_08_09_list)/2), 1):
#     H1N1_07_08_list.append(H1N1_08_09_list[i])
#     H1N1_07_08_filename_list.append(H1N1_08_09_filename_list[i])
#
# for i in range(0, int(len(H3N2_11_12_list)/2), 1):
#     H3N2_10_11_list.append(H3N2_11_12_list[i])
#     H3N2_10_11_filename_list.append(H3N2_11_12_filename_list[i])


H1N1_07_08_filename_list.sort()
H1N1_08_09_filename_list.sort()
H3N2_10_11_filename_list.sort()
H3N2_11_12_filename_list.sort()

# Organize the 4 subtype lists

# count = 0
# for i in range(int(len(H1N1_07_08_list)/2), len(H1N1_07_08_list), 1):
#     H1N1_07_08_list.insert(count, H1N1_07_08_list.pop(i))
#     H1N1_07_08_filename_list.insert(count, H1N1_07_08_filename_list.pop(i))
#     count += 2
count = 0
for i in range(int(len(H1N1_08_09_list)/2), len(H1N1_08_09_list), 1):
    H1N1_08_09_list.insert(count, H1N1_08_09_list.pop(i))
    H1N1_08_09_filename_list.insert(count, H1N1_08_09_filename_list.pop(i))
    count += 2
# count = 0
# for i in range(int(len(H3N2_10_11_list)/2), len(H3N2_10_11_list), 1):
#     H3N2_10_11_list.insert(count, H3N2_10_11_list.pop(i))
#     H3N2_10_11_filename_list.insert(count, H3N2_10_11_filename_list.pop(i))
#     count += 2
count = 0
for i in range(int(len(H3N2_11_12_list)/2), len(H3N2_11_12_list), 1):
    H3N2_11_12_list.insert(count, H3N2_11_12_list.pop(i))
    H3N2_11_12_filename_list.insert(count, H3N2_11_12_filename_list.pop(i))
    count += 2

count_subtype_list = [H5N1_list, H7N9_list, H1N1_list, H3N2_list,
                      H1N1_07_08_list, H1N1_08_09_list, H3N2_10_11_list, H3N2_11_12_list]
count_filename_list = [H5N1_filename_list, H7N9_filename_list, H1N1_filename_list, H3N2_filename_list,
                       H1N1_07_08_filename_list, H1N1_08_09_filename_list, H3N2_10_11_filename_list,
                       H3N2_11_12_filename_list]

output_file = open("Results/other_data.txt", "w")
print("Total number of unfiltered strains = " + str(counter_total_strains))
output_file.write("\nTotal number of unfiltered strains = " + str(counter_total_strains))
print("The number of strains containing ambiguous nucleotides = " + str(counter_ambiguous_strains))
output_file.write("\nThe number of strains containing ambiguous nucleotides = " + str(counter_ambiguous_strains))
print("The number of strains without a subtype = " + str(counter_no_subtype))
output_file.write("\nThe number of strains without a subtype = " + str(counter_no_subtype))
print("The number of duplicate strains = " + str(counter_duplicate_strains))
output_file.write("\nThe number of duplicate strains = " + str(counter_duplicate_strains))
# print("Number of strains after filtering = " + str(len(strain_list)))
print("Number of strains after filtering = " + str(counter_strains_post_filtering))
output_file.write("\nNumber of strains after filtering = " + str(counter_strains_post_filtering))
print()
output_file.write("\n")
print("The number of filtered strains per fasta file are shown below:")
output_file.write("\nThe number of filtered strains per fasta file are shown below:")

for i_1 in range(0, len(count_subtype_list), 1):
    for i_2 in range(0, len(count_subtype_list[i_1]), 1):
        temp_dict = count_subtype_list[i_1][i_2]
        temp_counter = len(temp_dict[count_filename_list[i_1][i_2]])
        print(str(count_filename_list[i_1][i_2]) + " = " + str(temp_counter))
        output_file.write("\n" + str(count_filename_list[i_1][i_2]) + " = " + str(temp_counter))

output_file.close()
H5N1_matrix_list = temp_strain.sum_codon_matrix_temp(H5N1_list, H5N1_filename_list)
H7N9_matrix_list = temp_strain.sum_codon_matrix_temp(H7N9_list, H7N9_filename_list)
H1N1_matrix_list = temp_strain.sum_codon_matrix_temp(H1N1_list, H1N1_filename_list)
H3N2_matrix_list = temp_strain.sum_codon_matrix_temp(H3N2_list, H3N2_filename_list)
H1N1_07_08_matrix_list = temp_strain.sum_codon_matrix_temp(H1N1_07_08_list, H1N1_07_08_filename_list)
H1N1_08_09_matrix_list = temp_strain.sum_codon_matrix_temp(H1N1_08_09_list, H1N1_08_09_filename_list)
H3N2_10_11_matrix_list = temp_strain.sum_codon_matrix_temp(H3N2_10_11_list, H3N2_10_11_filename_list)
H3N2_11_12_matrix_list = temp_strain.sum_codon_matrix_temp(H3N2_11_12_list, H3N2_11_12_filename_list)

H5N1_list = []
H7N9_list = []
H1N1_list = []
H3N2_list = []
H1N1_07_08_list = []
H1N1_08_09_list = []
H3N2_10_11_list = []
H3N2_11_12_list = []

final_list = []

for i in range(0, len(H5N1_matrix_list), 1):
    rscu = temp_strain.calc_rscu_value(H5N1_matrix_list[i])
    temp_list = [H5N1_filename_list[i]]
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    H5N1_list.append(temp_list)
for i in range(0, len(H7N9_matrix_list), 1):
    rscu = temp_strain.calc_rscu_value(H7N9_matrix_list[i])
    temp_list = [H7N9_filename_list[i]]
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    H7N9_list.append(temp_list)
for i in range(0, len(H1N1_matrix_list), 1):
    rscu = temp_strain.calc_rscu_value(H1N1_matrix_list[i])
    temp_list = [H1N1_filename_list[i]]
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    H1N1_list.append(temp_list)
for i in range(0, len(H3N2_matrix_list), 1):
    rscu = temp_strain.calc_rscu_value(H3N2_matrix_list[i])
    temp_list = [H3N2_filename_list[i]]
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    H3N2_list.append(temp_list)

# Seasonal Analysis
for i in range(0, len(H1N1_07_08_matrix_list), 1):
    rscu = temp_strain.calc_rscu_value(H1N1_07_08_matrix_list[i])
    temp_list = [H1N1_07_08_filename_list[i]]
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    H1N1_07_08_list.append(temp_list)
for i in range(0, len(H1N1_08_09_matrix_list), 1):
    rscu = temp_strain.calc_rscu_value(H1N1_08_09_matrix_list[i])
    temp_list = [H1N1_08_09_filename_list[i]]
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    H1N1_08_09_list.append(temp_list)
for i in range(0, len(H3N2_10_11_matrix_list), 1):
    rscu = temp_strain.calc_rscu_value(H3N2_10_11_matrix_list[i])
    temp_list = [H3N2_10_11_filename_list[i]]
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    H3N2_10_11_list.append(temp_list)
for i in range(0, len(H3N2_11_12_matrix_list), 1):
    rscu = temp_strain.calc_rscu_value(H3N2_11_12_matrix_list[i])
    temp_list = [H3N2_11_12_filename_list[i]]
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    H3N2_11_12_list.append(temp_list)

final_list.append(H5N1_list)
final_list.append(H7N9_list)
final_list.append(H1N1_list)
final_list.append(H3N2_list)
final_list.append(H1N1_07_08_list)
final_list.append(H1N1_08_09_list)
final_list.append(H3N2_10_11_list)
final_list.append(H3N2_11_12_list)

count = 0
for a in final_list:
    wb = Workbook()
    sheet_list = []
    for i in a:
        counter2 = 0
        for j in i:
            if isinstance(j, str):
                temp_sheet = wb.add_sheet(j)
            else:
                temp_sheet.write(counter2, 0, j[0])
                temp_sheet.write(counter2, 1, j[1])
                temp_sheet.write(counter2, 2, j[2])
            counter2 += 1
        sheet_list.append(temp_sheet)
    for i in sheet_list:
        i.write(0,0, "Codons")
        i.write(0,1, "Count")
        i.write(0,2, "RSCU Value")
    if count == 0:
        wb.save("Results/H5N1_Results.xls")
    elif count == 1:
        wb.save("Results/H7N9_Results.xls")
    elif count == 2:
        wb.save("Results/H1N1_Results.xls")
    elif count == 3:
        wb.save("Results/H3N2_Results.xls")
    elif count == 4:
        wb.save("Results/H1N1_07_08_Results.xls")
    elif count == 5:
        wb.save("Results/H1N1_08_09_Results.xls")
    elif count == 6:
        wb.save("Results/H3N2_10_11_Results.xls")
    elif count == 7:
        wb.save("Results/H3N2_11_12_Results.xls")
    count += 1





















