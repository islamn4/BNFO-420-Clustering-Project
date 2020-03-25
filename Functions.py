#Authour: Stephen Shea
#Created: 3/4/20
"""
This is a module that only contains functions and is intended for use in group 2's BNFO420 research project
"""


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


def H5N1_subtype_check(strain_header):
    """
        A function that determines if the virus strain in question is an H5N1 virus or a different subtype.

        :parameter
        strain_header : string
            Contains the fasta header of an Influenza A virus.
        :return:
        True if the virus strain is H5N1 strain.
        False if the virus strain is not an H5N1 strain.
        """

    if(strain_header.find("H5N1") != -1):
        return True
    else:
        return False


def subtype_filter(strain_header):
    """
        
    """

    if(strain_header.find("H5N1") != -1):
        return "H5N1"
    elif(strain_header.find("H9N2") != -1):
        return "H9N2"
    elif(strain_header.find("H1N1") != -1):
        return "H1N1"
    elif(strain_header.find("H3N2") != -1):
        return "H3N2"
    else:
        return "none"














































