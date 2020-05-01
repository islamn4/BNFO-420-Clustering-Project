class Strain(object):
    import NCBI_Fasta_Functions as ncbi
    """
        A class that creates an object for an Influenza-A virus strain.
    """

    EMPTY_CODON_DICT = {"UUU":0, "UUC":0, "UUA":0, "UUG":0, "CUU":0, "CUC":0, "CUA":0, "CUG":0, "AUU":0, "AUC":0,
                        "AUA":0, "AUG":0, "GUU":0, "GUC":0, "GUA":0, "GUG":0, "UCU":0, "UCC":0, "UCA":0, "UCG":0,
                        "CCU":0, "CCC":0, "CCA":0, "CCG":0, "ACU":0, "ACC":0, "ACA":0, "ACG":0, "GCU":0, "GCC":0,
                        "GCA":0, "GCG":0, "UAU":0, "UAC":0, "UAA":0, "UAG":0, "CAU":0, "CAC":0, "CAA":0, "CAG":0,
                        "AAU":0, "AAC":0, "AAA":0, "AAG":0, "GAU":0, "GAC":0, "GAA":0, "GAG":0, "UGU":0, "UGC":0,
                        "UGA":0, "UGG":0, "CGU":0, "CGC":0, "CGA":0, "CGG":0, "AGU":0, "AGC":0, "AGA":0, "AGG":0,
                        "GGU":0, "GGC":0, "GGA":0, "GGG":0}

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

    # Constructor method
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

    def __calc_rscu_summation__(cls, amino_acid, codon_matrix):
        """
        A class method that calculates the summation in the RSCU value equation.

        :param amino_acid: string
            Contains the single letter amino acid abrreviation for a specific codon
        :param codon_matrix:
        :return: int
            The denominator in the RSCU value equation.
        """
        summation = 0
        for i in cls.AA_DICT[amino_acid]:
            if isinstance(i, int):
                continue
            else:
                summation = summation + codon_matrix[i]
        return summation

    def calc_rscu_value(cls, codon_matrix = {}):
        """
            A class method that calculates the RSCU values of a given codon matrix
        :param codon_matrix:
        :return: dict

        """
        return_codon_matrix = {}
        temp_codon_dict = cls.CODON_DICT.copy()
        temp_codon_dict.pop("UAA")
        temp_codon_dict.pop("UAG")
        temp_codon_dict.pop("UGA")
        temp_codon_dict.pop("AUG")
        temp_codon_dict.pop("UGG")
        temp_aa_dict = cls.AA_DICT.copy()
        temp_aa_dict.pop("stop")
        temp_aa_dict.pop("W")
        temp_aa_dict.pop("M")
        for k,v in temp_codon_dict.items():
            amino_acid = v
            num_codons = temp_aa_dict[v][0]
            codon_freq = codon_matrix[k]
            summation = cls.__calc_rscu_summation__(amino_acid, codon_matrix)
            # if num_codons == 0:
            #     print("num_codons")
            # if summation == 0:
            #     print("summation")
            rscu = 0
            if num_codons != 0 and summation != 0:
                rscu = codon_freq / ((1 / num_codons) * summation)
            rscu = "{:.2f}".format(rscu)
            return_codon_matrix.update({k + "(" + v + ")":[codon_matrix[k], rscu]})
        return return_codon_matrix

    def sum_codon_matrix(self, matrix_list, total_strains):
        new_codon_matrix_list = []
        avian_ha_list = [["Avian_HA_H5N1"], ["Avian_HA_H7N9"]]
        avian_mp_list = [["Avian_MP_H5N1"], ["Avian_MP_H7N9"]]
        avian_na_list = [["Avian_NA_H5N1"], ["Avian_NA_H7N9"]]
        avian_pb2_list = [["Avian_PB2_H5N1"], ["Avian_PB2_H7N9"]]
        human_ha_list = [["Human_HA_H5N1"], ["Human_HA_H7N9"]]
        human_mp_list = [["Human_MP_H5N1"], ["Human_MP_H7N9"]]
        human_na_list = [["Human_NA_H5N1"], ["Human_NA_H7N9"]]
        human_pb2_list = [["Human_PB2_H5N1"], ["Human_PB2_H7N9"]]

        #Adds each strain objects to one of the eight lists above. Based on the host, gene, and subtype of the strain
        for i in matrix_list:
            # print(i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene())
            if ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H5N1_Avian_HA"):
                avian_ha_list[0].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H5N1_Avian_MP"):
                avian_mp_list[0].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H5N1_Avian_NA"):
                avian_na_list[0].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H5N1_Avian_PB2"):
                avian_pb2_list[0].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H5N1_Human_HA"):
                human_ha_list[0].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H5N1_Human_MP"):
                human_mp_list[0].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H5N1_Human_NA"):
                human_na_list[0].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H5N1_Human_PB2"):
                human_pb2_list[0].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H7N9_Avian_HA"):
                avian_ha_list[1].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H7N9_Avian_MP"):
                avian_mp_list[1].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H7N9_Avian_NA"):
                avian_na_list[1].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H7N9_Avian_PB2"):
                avian_pb2_list[1].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H7N9_Human_HA"):
                human_ha_list[1].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H7N9_Human_MP"):
                human_mp_list[1].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H7N9_Human_NA"):
                human_na_list[1].append(i)
                total_strains -= 1
            elif ((i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene()) == "H7N9_Human_PB2"):
                human_pb2_list[1].append(i)
                total_strains -= 1
            # else: print(i.get_subtype() + "_" + i.get_host() + "_" + i.get_gene() + " was not summed.")
        temp_strain_list = [avian_ha_list, avian_mp_list, avian_na_list, avian_pb2_list, human_ha_list, human_mp_list,
                            human_na_list, human_pb2_list]
        for i in temp_strain_list:
            counter = 0
            h5n1_codon_matrix = self.EMPTY_CODON_DICT.copy()
            h7n9_codon_matrix = self.EMPTY_CODON_DICT.copy()
            for h5n1 in i[0]:
                if isinstance(h5n1, str):
                    pass
                else:
                    for k in self.CODON_DICT.keys():
                        h5n1_codon_matrix[k] += h5n1.get_codon_matrix()[k]
                counter += 1
            counter = 0
            for h7n9 in i[1]:
                if isinstance(h7n9, str):
                    pass
                else:
                    for k in self.CODON_DICT.keys():
                        h7n9_codon_matrix[k] += h7n9.get_codon_matrix()[k]
                counter += 1

            h5n1_codon_matrix.pop("UAA")
            h7n9_codon_matrix.pop("UAA")

            h5n1_codon_matrix.pop("UAG")
            h7n9_codon_matrix.pop("UAG")

            h5n1_codon_matrix.pop("UGA")
            h7n9_codon_matrix.pop("UGA")

            h5n1_codon_matrix.pop("AUG")
            h7n9_codon_matrix.pop("AUG")

            h5n1_codon_matrix.pop("UGG")
            h7n9_codon_matrix.pop("UGG")

            i[0] = [i[0][0]]
            i[1] = [i[1][0]]
            i[0].append(h5n1_codon_matrix)
            i[1].append(h7n9_codon_matrix)
            new_codon_matrix_list.append(i[0])
            new_codon_matrix_list.append(i[1])

        # for i in new_codon_matrix_list:
        #     print()
        #     print(i)

        return new_codon_matrix_list

    def sum_codon_matrix_temp(self, strain_list, filename_list, total_strains = 0):
        matrix_list = []
        for i in range(0, len(strain_list), 1):
            temp_matrix = self.EMPTY_CODON_DICT.copy()
            for j in strain_list[i].get(filename_list[i]):
                for k in self.EMPTY_CODON_DICT.keys():
                    temp_matrix[k] += j.get_codon_matrix()[k]
            temp_matrix.pop("UAA")
            temp_matrix.pop("UAG")
            temp_matrix.pop("UGA")
            temp_matrix.pop("AUG")
            temp_matrix.pop("UGG")
            matrix_list.append(temp_matrix)
        # print()
        # print(matrix_list)
        return matrix_list

    def translate_dna_to_mrna(clss, seq):
        """

        :param seq:
        :return:
        """
        translated_seq = ""
        for char in seq:
            if char == "T":
                translated_seq = "A" + translated_seq
            elif char == "A":
                translated_seq = "U" + translated_seq
            elif char == "C":
                translated_seq = "G" + translated_seq
            elif char == "G":
                translated_seq = "C" + translated_seq
            else:
                print(char + " is an ambiguous nucleotide and cannot be translated.")

        return translated_seq

    def translate_mrna_to_dna(clss, seq):
        """

        :param seq:
        :return:
        """
        translated_seq = ""
        for char in seq:
            if char == "U":
                translated_seq = "A" + translated_seq
            elif char == "A":
                translated_seq = "T" + translated_seq
            elif char == "G":
                translated_seq = "C" + translated_seq
            elif char == "C":
                translated_seq = "G" + translated_seq
            else:
                print(char + " is an ambiguous nucleotide and cannot be translated.")
        return translated_seq


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

    def get_aa_dict(self):
        """

        :return:
        """
        return self.AA_DICT



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
        index = 0
        if self.ncbi.check_dna(seq):
            seq = seq[seq.find("ATG"):]
            for i in range(0, len(seq) + 1, 3):
                index += 3
                codon = seq[i:i + 3]
                if codon == "TAA":
                    seq = seq[:index]
                elif codon == "TAG":
                    seq = seq[:index]
                elif codon == "TGA":
                    seq = seq[:index]
                else: pass
            self.seq = seq
        else: print("The given sequence has ambiguous neucleotides or it is not a DNA sequence. "
                    "Please use a DNA sequence")

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
        seq = seq.replace("T", "U")
        codon_matrix = {}
        #
        for k in self.CODON_DICT.keys():
            codon_matrix.update({k: 0})
        # Iterates through the codons and counts the number of times each codon appears.
        for i in range(0, self.length + 1, 3):
            codon = seq[i:i+3]
            if(len(codon) == 3):
                if(codon_matrix.get(codon) == None):
                    print(codon + " does not exist in the codon matrix")
                else:
                    codon_matrix[codon] += 1
        self.codon_matrix = codon_matrix





def main():
    import csv
    test_file = open("test_files/HA_rscu_value_test_sequence.fasta", "r")
    seq = test_file.read()
    seq = seq[seq.find("\n"):].strip().replace("\n", "")

    test = Strain("HA" , "H5N1", "Avian", None, None, None, seq)

    # test.sum_codon_matrix([0])
    # print(test.calc_rscu_value(test.get_codon_matrix()))
    # print(len(test.get_seq()))
    temp_list = [["Codon", "Count", "RSCU"]]
    rscu = test.calc_rscu_value(test.get_codon_matrix())
    for k, v in rscu.items():
        temp_list.append([k, v[0], v[1]])
    with open('results.csv', 'w+') as f:
        writer = csv.writer(f)
        writer.writerows(temp_list)

if __name__ == "__main__":
    main()


