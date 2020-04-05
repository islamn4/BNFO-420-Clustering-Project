class Strain:
    import NCBI_Fasta_Functions as ncbi
    """
        A class that creates an object for an Influenza-A virus strain.
    """

    CODON_DICT = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "AUU":"I",
                  "AUC":"I", "AUA":"I", "AUG":"M", "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "UCU":"S", "UCC":"S",
                  "UCA":"S", "UCG":"S", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P", "ACU":"T", "ACC":"T", "ACA":"T",
                  "ACG":"T", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A", "UAU":"Y", "UAC":"Y", "UAA":"stop",
                  "UAG":"stop", "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                  "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "UGU":"C", "UGC":"C", "UGA":"stop", "UGG":"W", "CGU":"R",
                  "CGC":"R", "CGA":"R", "CGG":"R", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R", "GGU":"G", "GGC":"G",
                  "GGA":"G", "GGG":"G"}

    AA_DICT = {"F":[2,"UUU", "UUC"], "L":[6, "UUA","UUG","CUU","CUC","CUA","CUG"], "I":[3,"AUU","AUC","AUA"], "M":[1,"AUG"],
               "V":[4,"GUU","GUC","GUA","GUG"], "S":[6,"UCU","UCC","UCA","UCG","AGU","AGC"], "P":[4,"CCU","CCC","CCA","CCG"],
               "T":[4,"ACU","ACC","ACA","ACG"], "A":[4,"GCU","GCC","GCA","GCG"], "Y":[2,"UAU","UAC",], "H":[2,"CAU","CAC"],
               "Q":[2,"CAA","CAG"], "N":[2,"AAU","AAC"], "K":[2,"AAA","AAG"], "D":[2,"GAU","GAC"], "E":[2,"GAA","GAG"],
               "C":[2,"UGU","UGC"], "W":[1,"UGG"], "R":[6,"CGU","CGC","CGA","CGG","AGA","AGG"], "G":[4,"GGU","GGC","GGA","GGG"],
               "stop":[3,"UAA","UAG","UGA"]}

    #Constructor method
    def __init__(self, gene, subtype, host, country, year, accession_num, seq):
        self.set_gene(gene)
        self.set_subtype(subtype)
        self.set_host(host)
        self.set_country(country)
        self.set_year(year)
        self.set_accession_num(accession_num)
        self.set_seq(seq)
        self.set_length()
        self.set_codon_matrix()


    """
    # Constructor method
    def __init__(self, header, seq, host, filename):
        self.set_gene(gene)
        self.set_subtype(self.ncbi.parse_subtype(header))
        self.set_host(host)
        self.set_country(self.ncbi.parse_country(header))
        self.set_year(self.ncbi.parse_year(header))
        self.set_accession_num(self.ncbi.parse_id_num(header))
        self.set_seq(seq)
        self.set_length()
        self.set_codon_matrix()
    """
    def __str__(self):
        return self.accession_num

    def __calc_rscu_summation__(self,amino_acid, codon_matrix):
        """
        A function that calculates the summation in the RSCU value equation.

        :param num_codons: int
            Contains the number of codons that code for the provided amino acid.
        :param codon_freq: int
            Contains the number of times the provided amino acid is encoded by the provided codon
            for the sequence in question.
        :return: int
            The denominator in the RSCU value equation.
        """
        summation = 0
        for i in self.AA_DICT[amino_acid]:
            if isinstance(i, int):
                continue
            else:
                print(i)
                summation = summation + codon_matrix[i]
        return summation

    def calc_rscu_value(self, codon_matrix = {}):
        """

        :param codon_matrix:
        :return: dict

        """
        return_codon_matrix = {}
        for k,v in self.CODON_DICT.items():
            amino_acid = v
            num_codons = self.AA_DICT[v][0]
            codon_freq = codon_matrix[k]
            summation = self.__calc_rscu_summation__(amino_acid, codon_matrix)
            rscu = codon_freq / ((1 / num_codons) * summation)
            rscu = "{:.2f}".format(rscu)
            return_codon_matrix.update({k + "(" + v + ")":[codon_matrix[k], rscu]})
        return return_codon_matrix

    # def translate_dna_to_mrna(self, seq):
    #     """
    #
    #     :param seq:
    #     :return:
    #     """
    #     translated_seq = ""
    #     for char in seq:
    #         if char == "T":
    #             translated_seq = "A" + translated_seq
    #         elif char == "A":
    #             translated_seq = "U" + translated_seq
    #         elif char == "C":
    #             translated_seq = "G" + translated_seq
    #         elif char == "G":
    #             translated_seq = "C" + translated_seq
    #         else:
    #             print(char + " is an ambiguous nucleotide and cannot be translated.")
    #
    #     return translated_seq
    #
    # def translate_mrna_to_dna(self, seq):
    #     """
    #
    #     :param seq:
    #     :return:
    #     """
    #     translated_seq = ""
    #     for char in seq:
    #         if char == "U":
    #             translated_seq = "A" + translated_seq
    #         elif char == "A":
    #             translated_seq = "T" + translated_seq
    #         elif char == "G":
    #             translated_seq = "C" + translated_seq
    #         elif char == "C":
    #             translated_seq = "G" + translated_seq
    #         else:
    #             print(char + " is an ambiguous nucleotide and cannot be translated.")
    #     return translated_seq

    def get_gene(self):
        """

        :return:
        """
        return self.gene

    def get_subtype(self):
        """

        :return:
        """
        return self.subtype

    def get_host(self):
        """

        :return:
        """
        return self.host

    def get_country(self):
        """

        :return:
        """
        return self.country

    def get_year(self):
        """

        :return:
        """
        return self.year

    def get_accession_num(self):
        """

        :return:
        """
        return self.accession_num

    def get_seq(self):
        """

        :return: string
            The strain's DNA sequence.
        """
        return self.seq

    def get_length(self):
        """

        :return:
        """
        return self.length

    def get_codon_matrix(self):
        """
        Getter methond for instance variable codon_matrix
        :return: dictionary
            codon_matrix
        """
        return self.codon_matrix

    def get_codon_dict(self):
        """

        :return:
        """
        return self.CODON_DICT



    def set_gene(self, gene):
        """
        Setter method for instance variable gene.
        :param gene: string
            The gene that the nucleotide sequence codes for.
        :return:
            Nothing
        """
        self.gene = gene

    def set_host(self, host):
        """
        Sets the host species of the virus strain
        :param host: string
        :return:
            Nothing
        """
        self.host = host

    def set_subtype(self, subtype):
        """
        Setter method for instance variable subtype.
        :param subtype: string
            The subtype of the virus strain
        :return:
            Nothing
        """
        self.subtype = subtype

    def set_year(self, year):
        """
        Setter method for instance variable year.
        :param year: string
            The year the virus strain was isolated.
        :return:
            Nothing
        """
        self.year = year

    def set_accession_num(self, accession_num):
        """
        Setter method for instance variable accession_num.
        :param accession_num: string
            The virus strain's NCBI accession number.
        :return:
            Nothing
        """
        self.accession_num = accession_num

    def set_country(self, country):
        """
            Setter method for instance variable country.
        :param country: string
            The country where the virus strain was isolated.
        :return:
            Nothing
        """
        self.country = country

    def set_seq(self, seq):
        """
        Setter method for instance variable seq
        :param seq: string
            The virus strain's DNA and mRNA sequences
        :return:
            Nothing
        """
        seq = seq.upper().strip().replace("\n","")
        if self.ncbi.check_dna(seq):
            self.seq = seq
        else: print("The given sequence has ambiguous neucleotides or it is not a DNA sequence. "
                    "Please use a DNA or mRNA sequence")

    def set_length(self):
        """
        Setter method for instance variable length.
        :return:
            Nothing
        """
        self.length = int(len(self.get_seq()))

    def set_codon_matrix(self):
        """
            Setter method for instance variable codon_matrix
        :return: none
        """

        seq = self.get_seq()
        # seq = seq[seq.find("ATG"):]
        seq = seq.replace("T", "U")
        codon_matrix = {}
        #
        for k in self.CODON_DICT.keys():
            codon_matrix.update({k: 0})
        # Iterates through the codons
        for i in range(0, self.length + 1, 3):
            codon = seq[i:i+3]
            if(len(codon) == 3):
                if(codon_matrix.get(codon) == None):
                    print(codon + " does not exist in the codon matrix")
                    # pass
                else:
                    # print(codon_matrix.get(codon))
                    codon_matrix[codon] += 1
            # else: print(codon)
        self.codon_matrix = codon_matrix





def main():
    test_file = open("test_files/HA_rscu_value_test_sequence.fasta", "r")
    seq = test_file.read()
    seq = seq[seq.find("\n"):].strip().replace("\n", "")

    test = Strain(None, None, None, None, None, None, seq)


    print(test.calc_rscu_value(test.get_codon_matrix()))


if __name__ == "__main__":
    main()


