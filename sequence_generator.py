# 32 bp < input is needed

import re
import pandas as pd
import numpy as np
import xlsxwriter
import random
import openpyxl

SEQUENCE_LENGTH_FOR_VALIDATION = 32
RANDOM_SAMPLING_CONSTANT = 60
GUIDE_LENGTH = 20
NUMBER_OF_FULLY_MATCHED_TARGET = 60
nt = ["A", "T", "G", "C"]
NT = 4

# input validation
def seq_validator(data):
    if data[0] == "#":
        return None
    m = re.findall(r"^[A|a|T|t|C|c|G|g]+$", data)
    return m[0] if m else None


# input pre-processing
def input_processor(file):
    """_summary_
    dictionary element of {gene name : 32 bp sequence (context + PAM + protospacer)

    template

    Args:
        file (_type_): _description_

    Returns:
        _type_: _description_
    """
    df = pd.read_excel(file, engine="openpyxl", index_col=0)
    seq_collection = df.set_index(["#Gene Name"])[
        "#Sequence with context (4 + 4 + 20 + 4)"
    ].to_dict()
    df = None

    # dictionary pre-processing
    keys_to_deleted = []
    for item in seq_collection:
        if seq_validator(seq_collection[item]) is None:
            keys_to_deleted.append(item)

    for key in keys_to_deleted:
        seq_collection.pop(key)

    return seq_collection


class Mutation:
    # mutation dictionary
    ALL_MUTATION = {
        "A": ["g", "c", "t"],
        "G": ["a", "c", "t"],
        "C": ["a", "g", "t"],
        "T": ["a", "g", "c"],
    }
    TRANSITION = {"A": ["g"], "G": ["a"], "C": ["t"], "T": ["c"]}
    TRANSVERSION = {"A": ["c", "t"], "G": ["c", "t"], "C": ["a", "g"], "T": ["a", "g"]}
    TRANSVERSION_SIM_GUIDE = {"A": "t", "G": "c", "C": "g", "T": "a"}

    def all_mutation(self, seed):
        """_summary_

        Args:
            seed (_type_): _description_

        Returns:
            single mutation only
        """
        df = pd.DataFrame(
            columns=[
                "Original sequence",
                "mutant sequence",
                "gene name",
                "mutation type",
            ]
        )

        for gene_name in seed:
            try:

                # TODO: move to input processor
                # input sequence length qualification
                if len(seed[gene_name]) < SEQUENCE_LENGTH_FOR_VALIDATION:
                    print("Sequence information is invalid. The operation is aborted.")
                    return -1

                # input seed pre-processing (partitioning)
                lead_4bp = seed[gene_name][:4]
                original_pam = seed[gene_name][4:8]
                protospacer = seed[gene_name][8:28]
                trailing_4bp = seed[gene_name][28:]

                # making a 1-bp mismatch collection for a single guide sequence
                mt_protospacer_collection = list()
                for i in range(len(protospacer)):
                    for occurrence in self.ALL_MUTATION[protospacer[i]]:
                        mt_protospacer_collection.append(
                            f"{protospacer[:i]}{occurrence}{protospacer[i+1:]}"
                        )

                for entry in mt_protospacer_collection:
                    df = df.append(
                        {
                            "Original sequence": f"{seed[gene_name]}",
                            "mutant sequence": f"{lead_4bp}{original_pam}{entry}{trailing_4bp}",
                            "gene name": f"{gene_name}",
                            "mutation type": "ALL_POINT_MUTATION",
                        },
                        ignore_index=True,
                    )
            except Exception as e:
                print(e)
                print("PAM library generation aborted")
                raise

        return df

    def transition(self, seed, n=0):
        # # = n mutation functionality not implemented yet
        """_summary_

        Args:
            seed (_type_): _description_

        Returns:
            _type_: _description_
        """
        if n == 0:
            print("invalid input. Operation aborted")
            return -1

        df = pd.DataFrame(
            columns=[
                "Original sequence",
                "mutant sequence",
                "gene name",
                "mutation type",
            ]
        )

        for gene_name in seed:
            try:

                # TODO: move to input processor
                # input sequence length qualification
                if len(seed[gene_name]) < SEQUENCE_LENGTH_FOR_VALIDATION:
                    print("Sequence information is invalid. The operation is aborted.")
                    return -1

                # input seed pre-processing (partitioning)
                lead_4bp = seed[gene_name][:4]
                original_pam = seed[gene_name][4:8]
                protospacer = seed[gene_name][8:28]
                trailing_4bp = seed[gene_name][28:]

                # making a 1-bp mismatch collection for a single guide sequence
                mt_protospacer_collection = list()
                for i in range(len(protospacer)):
                    for occurrence in self.TRANSITION[protospacer[i]]:
                        mt_protospacer_collection.append(
                            f"{protospacer[:i]}{occurrence}{protospacer[i+1:]}"
                        )

                for entry in mt_protospacer_collection:
                    df = df.append(
                        {
                            "Original sequence": f"{seed[gene_name]}",
                            "mutant sequence": f"{lead_4bp}{original_pam}{entry}{trailing_4bp}",
                            "gene name": f"{gene_name}",
                            "mutation type": "TRANSITION",
                        },
                        ignore_index=True,
                    )
            except Exception as e:
                print(e)
                print("PAM library generation aborted")
                raise

        return df

    def transversion(self, seed, n=0):
        """_summary_

        Args:
            seed (_type_): _description_

        Returns:
            _type_: _description_
        """

        if n == 0:
            print("invalid input. Operation aborted")
            return -1
        df = pd.DataFrame(
            columns=[
                "Original sequence",
                "mutant sequence",
                "gene name",
                "mutation type",
            ]
        )

        for gene_name in seed:
            try:

                # TODO: move to input processor
                # input sequence length qualification
                if len(seed[gene_name]) < SEQUENCE_LENGTH_FOR_VALIDATION:
                    print("Sequence information is invalid. The operation is aborted.")
                    return -1

                # input seed pre-processing (partitioning)
                lead_4bp = seed[gene_name][:4]
                original_pam = seed[gene_name][4:8]
                protospacer = seed[gene_name][8:28]
                trailing_4bp = seed[gene_name][28:]

                # generate all possible combination of n-bp mismatches
                mt_protospacer_collection = list()

                # generating serial mutants first, when multiple mutations needed
                if n > 1:
                    for pos in range(GUIDE_LENGTH - (n - 1)):
                        if n == 2:
                            i, j = pos, pos + 1
                            temp_mt = f"{protospacer[:i]}{self.TRANSVERSION_SIM_GUIDE[protospacer[i]]}{protospacer[i+1:j]}{self.TRANSVERSION_SIM_GUIDE[protospacer[j]]}{protospacer[j+1:]}"
                            mt_protospacer_collection.append(temp_mt)

                        elif n == 3:
                            i, j, k = pos, pos + 1, pos + 2
                            temp_mt = f"{protospacer[:i]}{self.TRANSVERSION_SIM_GUIDE[protospacer[i]]}{protospacer[i+1:j]}{self.TRANSVERSION_SIM_GUIDE[protospacer[j]]}{protospacer[j+1:k]}{self.TRANSVERSION_SIM_GUIDE[protospacer[k]]}{protospacer[k+1:]}"
                            mt_protospacer_collection.append(temp_mt)
                        elif n == 4:
                            i, j, k, l = pos, pos + 1, pos + 2, pos + 3
                            temp_mt = f"{protospacer[:i]}{self.TRANSVERSION_SIM_GUIDE[protospacer[i]]}{protospacer[i+1:j]}{self.TRANSVERSION_SIM_GUIDE[protospacer[j]]}{protospacer[j+1:k]}{self.TRANSVERSION_SIM_GUIDE[protospacer[k]]}{protospacer[k+1:l]}{self.TRANSVERSION_SIM_GUIDE[protospacer[l]]}{protospacer[l:]}"
                            mt_protospacer_collection.append(temp_mt)

                # generating random mutants
                while len(mt_protospacer_collection) < RANDOM_SAMPLING_CONSTANT * (
                    2 * n - 1
                ):

                    # making a 1-bp mismatch collection for a single guide sequence
                    if n == 1:
                        for i in range(len(protospacer)):
                            for occurrence in self.TRANSVERSION_SIM_GUIDE[
                                protospacer[i]
                            ]:
                                mt_protospacer_collection.append(
                                    f"{protospacer[:i]}{occurrence}{protospacer[i+1:]}"
                                )
                        break
                    elif n == 2:
                        indices = random.sample(range(0, GUIDE_LENGTH), n)
                        # serial assignment
                        indices.sort()
                        i, j = indices
                        temp_mt = f"{protospacer[:i]}{self.TRANSVERSION_SIM_GUIDE[protospacer[i]]}{protospacer[i+1:j]}{self.TRANSVERSION_SIM_GUIDE[protospacer[j]]}{protospacer[j+1:]}"
                        if temp_mt not in mt_protospacer_collection:
                            mt_protospacer_collection.append(temp_mt)

                    elif n == 3:
                        indices = random.sample(range(0, GUIDE_LENGTH), n)
                        # serial assignment
                        indices.sort()
                        i, j, k = indices
                        temp_mt = f"{protospacer[:i]}{self.TRANSVERSION_SIM_GUIDE[protospacer[i]]}{protospacer[i+1:j]}{self.TRANSVERSION_SIM_GUIDE[protospacer[j]]}{protospacer[j+1:k]}{self.TRANSVERSION_SIM_GUIDE[protospacer[k]]}{protospacer[k+1:]}"
                        if temp_mt not in mt_protospacer_collection:
                            mt_protospacer_collection.append(temp_mt)
                    elif n == 4:
                        indices = random.sample(range(0, GUIDE_LENGTH), n)
                        # serial assignment
                        indices.sort()
                        i, j, k, l = indices
                        temp_mt = f"{protospacer[:i]}{self.TRANSVERSION_SIM_GUIDE[protospacer[i]]}{protospacer[i+1:j]}{self.TRANSVERSION_SIM_GUIDE[protospacer[j]]}{protospacer[j+1:k]}{self.TRANSVERSION_SIM_GUIDE[protospacer[k]]}{protospacer[k+1:l]}{self.TRANSVERSION_SIM_GUIDE[protospacer[l]]}{protospacer[l:]}"
                        if temp_mt not in mt_protospacer_collection:
                            mt_protospacer_collection.append(temp_mt)
                    # display current generation status
                    print(
                        f"{len(mt_protospacer_collection)}-th unique mutation has been generated."
                    )

                # save to Pandas dataframe
                for entry in mt_protospacer_collection:
                    df = df.append(
                        {
                            "Original sequence": f"{seed[gene_name]}",
                            "mutant sequence": f"{lead_4bp}{original_pam}{entry}{trailing_4bp}",
                            "gene name": f"{gene_name}",
                            "mutation type": "TRANSVERSION",
                        },
                        ignore_index=True,
                    )
            except Exception as e:
                print(e)
                print("PAM library generation aborted")
                raise

        return df

    # Mismatched target (1bp)
    def mt_1(self, seed):
        """mismatched for a single position ("mt" refers to the term, mutant or mismatch)

        Args:
            seed (dict): gene name : guide sequences (4+4+20+4 bp, 5' to 3')
        """
        df = self.all_mutation(seed)
        # Output as an excel file
        writer = pd.ExcelWriter("Mismatched target (1bp).xlsx", engine="xlsxwriter")
        df.to_excel(writer, sheet_name="1 bp MM")
        writer.save()

    # Mismatched target (2bp)
    def mt_2(self, seed):
        """_summary_
            transversion
        Args:
            seed (_type_): _description_
        """

        df = self.transversion(seed, 2)
        # Output as an excel file
        writer = pd.ExcelWriter("Mismatched target (2bp).xlsx", engine="xlsxwriter")
        df.to_excel(writer, sheet_name="2 bp MM")
        writer.save()

    # Mismatched target (3bp)
    def mt_3(self, seed):
        """_summary_
            transversion
        Args:
            seed (_type_): _description_
        """

        df = self.transversion(seed, 3)
        # Output as an excel file
        writer = pd.ExcelWriter("Mismatched target (3bp).xlsx", engine="xlsxwriter")
        df.to_excel(writer, sheet_name="3 bp MM")
        writer.save()

    # Mismatched target (4bp)
    def mt_4(self, seed):
        """_summary_
            transversion
        Args:
            seed (_type_): _description_
        """

        df = self.transversion(seed, 4)
        # Output as an excel file
        writer = pd.ExcelWriter("Mismatched target (4bp).xlsx", engine="xlsxwriter")
        df.to_excel(writer, sheet_name="4 bp MM")
        writer.save()


# PAM variants
def PAM(seed, length=4):
    """Generate all possible oligonucleotide sequence with all PAM sequence for Cas12f1 (NNNN)

    making excel spreadsheet

    Input sequence length (format) should be considered carefully.
    4 bp(trailing) + 4 bp (PAM) + 20 bp (guide) + 4 bp (trailing) + a (+alpha) will be rendered.

    Args:
        seed (dict): gene name : guide sequences (4+4+20+4 bp, 5' to 3')
        length (int, optional): _description_. Defaults to 4.


    """
    if not isinstance(seed, dict):
        seed = dict(seed)

    pam256 = list()
    # generate all possible PAM combinations
    # O(n^4) : not suitable for long PAM sequence
    for i in range(NT):
        for j in range(NT):
            for k in range(NT):
                for l in range(NT):
                    pam256.append(f"{nt[i]}{nt[j]}{nt[k]}{nt[l]}")

    df_PAM = pd.DataFrame(
        columns=[
            "PAM",
            "Gene Name",
            "Seed guide sequence",
            "PAM + guide",
            "32 bp synthetic target and target context sequence (4 bp + PAM + 20 bp protospacer + 4 bp)",
        ]
    )
    for p in pam256:
        for gene_name in seed:
            try:

                # TODO: move to input processor
                # input sequence length qualification
                if len(seed[gene_name]) < SEQUENCE_LENGTH_FOR_VALIDATION:
                    print("Sequence information is invalid. The operation is aborted.")
                    return -1

                # input seed pre-processing
                lead_4bp = seed[gene_name][:4]
                original_pam = seed[gene_name][4:8]
                protospacer = seed[gene_name][8:28]
                trailing_4bp = seed[gene_name][28:]

                df_PAM = df_PAM.append(
                    {
                        "PAM": f"{p}",
                        "Gene Name": f"{gene_name}",
                        "Seed guide sequence": f"{protospacer}",
                        "PAM + guide": f"{p}{protospacer}",
                        "32 bp synthetic target and target context sequence (4 bp + PAM + 20 bp protospacer + 4 bp)": f"{lead_4bp}{p}{protospacer}{trailing_4bp}",
                    },
                    ignore_index=True,
                )
            except Exception as e:
                print(e)
                print("PAM library generation aborted")
                raise

    # Output as an excel file
    writer = pd.ExcelWriter("PAM Library.xlsx", engine="xlsxwriter")
    df_PAM.to_excel(writer, sheet_name="PAM Library")
    writer.save()


def guide_random_generator(iter):
    """_summary_"""
    guide_collection = pd.DataFrame(
        columns=[
            "20 bp Randomly generated fully matched targets",
        ]
    )

    # Randomly generating unique guide sequences
    for i in range(int(iter)):
        string_temp = "TTTR"

        for j in range(GUIDE_LENGTH - 4):
            string_temp += nt[random.sample(range(0, len(nt)), 1)[0]]

        # TODO: is there any better way
        # no duplicates
        if (
            guide_collection["20 bp Randomly generated fully matched targets"]
            .str.contains(string_temp)
            .any()
        ):
            continue
        guide_collection = guide_collection.append(
            {"20 bp Randomly generated fully matched targets": f"{string_temp}"},
            ignore_index=True,
        )
        print(f"{i}-th random sequence has been generated.")

    # Output as an excel file
    writer = pd.ExcelWriter(
        "20 bp Randomly generated fully matched targets.Xlsx", engine="xlsxwriter"
    )
    guide_collection.to_excel(writer, sheet_name="random generations")
    writer.save()


# Fully-matched target
# Guide sequence length


# Start of the program
"""
if __name__ == "__main__":
    # Takes input of seed guides
    # fixed name of inputs (guide.txt)

    # Which actions will be executed is determined by the user
    input("Choose the action ...")

    with open("guide.txt",'r') as file:
        for line in file:
            # do necessary processing line by line
            print(line) 


    PAM(obj)
"""

# input format
# sequence and name of that are separated by a white space


# debugging area #############################################

# PAM({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"})


# print(mt.all_mutation({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"}))
# print(mt.transition({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"}))
# print(mt.transversion({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"}, 2))
# print(mt.transversion({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"}, 3))
# print(mt.transversion({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"}, 4))
# mt.mt_1({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"})
# mt.mt_2({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"})
# mt.mt_3({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"})
# mt.mt_4({"LIN7B": "TCTATTTAACTGCCTAGTCGAACTCAATATCTCC"})
# guide_random_generator(10)

src = "input_template.xlsx"

# hard-coded operations
mt = Mutation()
PAM(input_processor(src))
mt.mt_1(input_processor(src))
mt.mt_2(input_processor(src))
mt.mt_3(input_processor(src))
mt.mt_4(input_processor(src))
guide_random_generator(NUMBER_OF_FULLY_MATCHED_TARGET)
