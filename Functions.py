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











































