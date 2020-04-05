#Author: Stephen Shea
#Created: 4/2/20
"""
    The functions return specific parts of a single NCBI fasta header.
"""
import re

def check_ambiguous_bases(seq):
    """
    A function that looks for ambiguous nucleotides.
    The ambiguous nucleotides include, without the quotes, the following:
        "N", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V".

    :parameter
    seq : string
        Contains the nucleotide sequence to be searched.
    :return:
    True if at least one ambiguous nucleotide is found.
    False if no ambiguous nucleotides are found.
    """

    seq = seq.upper()
    if re.search("[^ATGCU\s]", seq) == None:
        return False
    else:
        return True

def check_ambiguous_aa(seq):
    """
    A function that looks for ambiguous amino acids.
    The ambiguous amino acids include, without the quotes, the following: "B", "Z", "X", "J", "O", "U", ".", "*"

    :parameter
    seq : string
        Contains the amino acid sequence to be searched
    :return:
    True if if at least one ambiguous amino acids are found.
    False if no ambiguous amino acid is found.
    """

    seq = seq.upper()

    if(seq.find("B") != -1):
        return True
    elif(seq.find("Z") != -1):
        return True
    elif(seq.find("X") != -1):
        return True
    elif(seq.find("J") != -1):
        return True
    elif(seq.find("O") != -1):
        return True
    elif (seq.find("U") != -1):
        return True
    elif(seq.find(".") != -1):
        return True
    elif(seq.find("*") != -1):
        return True
    else:
        return False

def check_mrna(seq):
    """
    Function that checks if the given sequence is an mRNA sequence.
    :param seq: string
        Contains the sequence to be searched.
    :return:
        True if the sequence is an mRNA sequence.
        False if the sequence is not an mRNA sequence.
    """
    seq = seq.upper()
    if re.search("[^AGCU\s]", seq) == None:
        return True
    else:
        return False

def check_dna(seq):
    """
        Function that checks if the given sequence is a DNA sequence.
        :param seq: string
            Contains the sequence to be searched.
        :return:
            True if the sequence is a DNA sequence.
            False if the sequence is not an DNA sequence.
        """
    seq = seq.upper()
    if re.search("[^AGCT\s]", seq) == None:
        return True
    else:
        return False

def parse_id_num(header):
    """
    Returns the accession number from the header
    :param header: string
        NCBI fasta header
    :return: string
        NCBI accession number for the sequence
    """
    if(header.find(">") == 0):
        header = header[1:]
    id_num = header[:header.find(" ") + 1]

    return id_num

def parse_subtype(header):
    """

    :param header:
    :return:
    """
    subtype = re.search("[hH]\d[nN]\d", header)
    if subtype != None:
        return subtype.group(0)
    else: return subtype

# def parse_host(header, filename = False):
#     if filename:
#
#     else:

def parse_year(header):
    """

    :param header:
    :return:
    """
    year = re.search("[\/][20]\d+[(][hH]", header)
    if year != None:
        year = year.group(0)
        year = year[1:year.find("(")]
        if year[0] == "0":
            year = "20" + year
        # print(year)
    return year

#Needs more testing
def parse_country(header):
    """

    :param header:
    :return:
    """
    country = re.search("[A][\/].+[(][H]\d[N]\d[)]", header)
    if country != None:
        country = country.group(0)
        temp_list = country.split("/")
        country = temp_list.pop(len(temp_list) - 3)
        # print(country)
    return(country)

#Does not work yet
def get_gene(header):
    gene = "none"
    ha = re.search("[(][H][A][)]", header)
    # if(ha == None):
    #     ha = re.match(" segment 4 ", header)
    if(ha != None):
        gene = ha.group(0)[1:3]
    elif(re.match(" segment 4 ", header) != None):
        gene = "HA"
    elif((re.match("M1") != None) and (re.match("M2") != None)):
        gene = "M1_M2"

































