#Author: Stephen Shea
#Created: 4/12/20
"""
    The functions return specific parts of a single Flu DB fasta header.
"""
import re


def parse_id_num(header):
    """
    Returns the accession number from the header
    :param header: string
        Flu_DB fasta header
    :return: string
        Flu_DB accession number for the sequence
    """
    if header.find(">") == 0:
        header = header[1:]
    id_num = header[header.find(":")+1:header.find("|")]
    return id_num

def parse_year(header):
    """

    :param header:
    :return:
    """
    year = re.search("[\/][1290]\d+[|][pP]", header)
    # if len(year) <= 3:
    #     year = re.search("[\/][120]\d+[(][hH]", header)

    if year != None:
        year = year.group(0)
        year = year[1:year.find("|")].strip()
        if len(year) <= 3:
            # print()
            # print(header)
            year = re.search("[\/][120]\d+[(][hH]", header)
            year = year.group(0)
            year = year[1:year.find("(")]
        if year[0] == "0":
            year = "20" + year
        # print(year)
    return year

def parse_gene(header):
    """

    :param header:
    :return:
    """
    if header.find(">") == 0:
        header = header[1:]
    gene = header[header.find(":", header.find("Gene Symbol:"))+1:header.find("|", header.find("Gene Symbol:"))]
    gene = gene.strip().upper()
    return gene




# def check_duplicates(temp_strain, strain_list, filename):
#     """
#
#     :param strain:
#     :param strain_list:
#     :return:
#     """
#
#     for strain in strain_list[len(strain_list) - 1].get(str(filename)):
#         if temp_strain.get_accession_num() == strain.get_accession_num():
#             counter_duplicate_strains += 1
#             duplicate = True
#             break
#         elif temp_strain.get_seq() == strain.get_seq():
#             counter_duplicate_strains += 1
#             duplicate = True
#             break
#     if not duplicate:
