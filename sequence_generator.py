# 32 bp < input is needed
import itertools
import random
import re

import pandas as pd

# pre-defined values
SEQUENCE_LENGTH_FOR_VALIDATION = 32
RANDOM_SAMPLING_CONSTANT = 60
parameter = 20
NUMBER_OF_FULLY_MATCHED_TARGET = 60
PAM_FOR_RANDOM_GENERTION = ["TTTA", "TTTG"]
nt = ["A", "T", "G", "C"]
NT = 4


class guide_parameter():

    # user input for setting the guide parameters
    def __init__(self):
        print("Write proper information on your dataset")
        self.len_leading_seq = int(input("leading sequence length: "))
        self.len_pam = int(input("PAM length: "))
        self.len_guide_seq = int(input("guide sequence length: "))
        self.len_trailing_seq = int(input("trailing sequence length: "))
        self.number_of_maximum_mismatches = int(input("maximum number of mismatches appeared in the sequence: "))

        # debug
        # self.len_leading_seq = 150
        # self.len_pam = 4
        # self.len_guide_seq = 20
        # self.len_trailing_seq = 126
        # self.number_of_maximum_mismatches = 4


# input validation
def seq_validator(data):
    if data[0] == "#":
        return None
    m = re.findall(r"^[AaTtCcGg]+$", data)
    return m[0] if m else None


# input pre-processing
def input_processor(file):
    """_summary_
    Read excel file, and transfer it to Pandas DataFrame

    template

    Args:
        file (_type_): _description_

    Returns:
        _type_: _description_
    """
    df = pd.read_excel(file, engine="openpyxl")
    df.astype({"#Indel frequency": "int32"})

    # TODO: how to determine whether the input is valid or not
    # input sequence length validation

    return df


def sequence_partitioner(row, parameter: guide_parameter):
    leading_sequence = row["#Sequence with context (150 + 4 + 20 + 126)"][:parameter.len_leading_seq]
    original_pam = row["#Sequence with context (150 + 4 + 20 + 126)"][
                   parameter.len_leading_seq:parameter.len_leading_seq + parameter.len_pam]
    protospacer = row["#Sequence with context (150 + 4 + 20 + 126)"][
                  parameter.len_leading_seq + parameter.len_pam:parameter.len_leading_seq + parameter.len_pam + parameter.len_guide_seq]
    trailing_sequence = row["#Sequence with context (150 + 4 + 20 + 126)"][
                        parameter.len_leading_seq + parameter.len_pam + parameter.len_guide_seq:]

    return leading_sequence, original_pam, protospacer, trailing_sequence


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

    def all_mutation(self, df_input_sequence, parameter):
        """_summary_

        Args:
            parameter:
            df_input_sequence (_type_): _description_

        Returns:
            single mutation only
        """
        df = pd.DataFrame(
            columns=[
                "Original sequence",
                "mutant guide sequence",
                "Generated sequence with gDNA context",
                "gene name",
                "Known indel frequency from the source",
                "mutation information",
                "mutation type",
            ]
        )
        skipped = int()
        for index, row in df_input_sequence.iterrows():
            try:
                if pd.isna(row["#Sequence with context (150 + 4 + 20 + 126)"]):
                    skipped += 1
                    continue
                # input seed pre-processing (partitioning)
                leading_sequence, original_pam, protospacer, trailing_sequence = sequence_partitioner(row, parameter)

                # making a 1-bp mismatch collection for a single guide sequence
                mt_protospacer_collection = list()
                mt_info_collection = list()
                for i in range(len(protospacer)):
                    for occurrence in self.ALL_MUTATION[protospacer[i]]:
                        mt_protospacer_collection.append(
                            f"{protospacer[:i]}{occurrence}{protospacer[i + 1:]}"
                        )
                        # (loc, a -> b)
                        mt_info_collection.append((i, f"{protospacer[i]} -> {occurrence}"))

                for entry, supple_info in zip(mt_protospacer_collection, mt_info_collection):
                    df = df.append(
                        {
                            "Original sequence": f"{row['#Sequence with context (150 + 4 + 20 + 126)']}",
                            "mutant guide sequence": f"{entry}",
                            "Generated sequence with gDNA context": f"{leading_sequence}{original_pam}{entry}{trailing_sequence}",
                            "mutation information": f"location on the guide: {supple_info[0]}, {supple_info[1]}",
                            "gene name": f"{row['#Gene Name']}",
                            "Known indel frequency from the source": f"{row['#Indel frequency']}",
                            "mutation type": "ALL_POINT_MUTATION",
                        },
                        ignore_index=True,
                    )
            except Exception as e:
                print(e)
                print("PAM library generation aborted")
                raise

        print(f"{skipped} occurrences has been skipped")
        return df

    def transversion(self, df_input_sequence, parameter: guide_parameter, num_mismatches=2):
        """_summary_

        Args:
            num_mismatches:
            df_input_sequence (_type_): _description_

        Returns:
            _type_: _description_
        """
        skipped = int()

        df = pd.DataFrame(
            columns=[
                "Original sequence",
                "mutant guide sequence",
                "Generated sequence with gDNA context",
                "mutation information",
                "gene name",
                "mutation type",
                "Known indel frequency from the source",
            ]
        )
        for index, row in df_input_sequence.iterrows():
            try:
                if pd.isna(row["#Sequence with context (150 + 4 + 20 + 126)"]):
                    skipped += 1
                    continue

                # input seed pre-processing (partitioning)
                leading_sequence, original_pam, protospacer, trailing_sequence = sequence_partitioner(row, parameter)

                # generate all possible combination of n-bp mismatches
                mt_protospacer_collection = list()
                mt_info_collection = list()

                tokens = list()
                # generating serial mutants first, when multiple mutations needed
                for pos in range(parameter.len_guide_seq - num_mismatches + 1):
                    # generating indices
                    indices = [pos + i for i in range(num_mismatches)]

                    temp_mut_info_list = list()

                    # stringify list elements
                    mutant_string = str()
                    mt_info_string = str()

                    temp_list = list(protospacer)
                    for idx in indices:
                        temp_list[idx] = self.TRANSVERSION_SIM_GUIDE[protospacer[idx]]
                        # (loc, a -> b)
                        temp_mut_info_list.append((
                            f"{idx}",
                            f"{protospacer[idx]} -> {self.TRANSVERSION_SIM_GUIDE[protospacer[idx]]}"))

                    mutant_string = "".join(temp_list)
                    # (loc, a -> b)
                    for entry in temp_mut_info_list:
                        mt_info_string += f"{entry[0]}, {entry[1]}\n"

                    # freeing memory
                    temp_mut_info_list = None

                    # store processed string in a list
                    mt_protospacer_collection.append(mutant_string)
                    mt_info_collection.append(mt_info_string)

                # generating random mutants
                while len(mt_protospacer_collection) < RANDOM_SAMPLING_CONSTANT * (
                        2 * num_mismatches - 1
                ):
                    # making a 1-bp mismatch collection for a single guide sequence, default case?
                    if parameter.number_of_maximum_mismatches == 1:
                        for i in range(len(protospacer)):
                            for occurrence in self.TRANSVERSION_SIM_GUIDE[
                                protospacer[i]
                            ]:
                                mt_protospacer_collection.append(
                                    f"{protospacer[:i]}{occurrence}{protospacer[i + 1:]}"
                                )
                        break
                    else:
                        # randomly generated mutation locations
                        indices = random.sample(range(0, parameter.len_guide_seq),
                                                num_mismatches)
                        # serial assignment
                        indices.sort()

                        temp_mut_info_list = list()

                        # stringify list elements
                        mutant_string = str()
                        mt_info_string = str()

                        temp_list = list(protospacer)
                        for idx in indices:
                            temp_list[idx] = self.TRANSVERSION_SIM_GUIDE[protospacer[idx]]
                            # (loc, a -> b)
                            temp_mut_info_list.append((
                                f"{idx}",
                                f"{protospacer[idx]} -> {self.TRANSVERSION_SIM_GUIDE[protospacer[idx]]}"))

                        mutant_string = "".join(temp_list)
                        # (loc, a -> b)
                        for entry in temp_mut_info_list:
                            mt_info_string += f"{entry[0]}, {entry[1]}\n"

                        # freeing memory
                        temp_mut_info_list = None

                        # store processed string in a list only if it is unique
                        if mutant_string not in mt_protospacer_collection:
                            mt_protospacer_collection.append(mutant_string)
                            mt_info_collection.append(mt_info_string)

                    # display current generation status
                    print(
                        f"{len(mt_protospacer_collection)}-th unique mutation has been generated."
                    )

                # save to Pandas dataframe
                for entry, supp_info in zip(mt_protospacer_collection, mt_info_collection):
                    df = df.append(
                        {
                            "Original sequence": f"{row['#Sequence with context (150 + 4 + 20 + 126)']}",
                            "mutant guide sequence": entry,
                            "Generated sequence with gDNA context": f"{leading_sequence}{original_pam}{entry}{trailing_sequence}",
                            "mutation information": supp_info,
                            "gene name": f"{row['#Gene Name']}",
                            "mutation type": "TRANSVERSION",
                            "Known indel frequency from the source": f"{row['#Indel frequency']}",

                        },
                        ignore_index=True,
                    )
            except Exception as e:
                print(e)
                print("PAM library generation aborted")
                raise
        print(f"{skipped} ocurrences has been skipped")
        return df

    '''
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
                "mutant guide sequence",
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
                            f"{protospacer[:i]}{occurrence}{protospacer[i + 1:]}"
                        )

                for entry in mt_protospacer_collection:
                    df = df.append(
                        {
                            "Original sequence": f"{seed[gene_name]}",
                            "mutant guide sequence": f"{lead_4bp}{original_pam}{entry}{trailing_4bp}",
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
    '''

    # Mismatched target (1bp)
    def mt(self, df_input_sequence, parameter: guide_parameter):
        """mismatched for a single position ("mt" refers to the term, mutant or mismatch)
            for multiple mismatches, the program only considers transversion mutations
        Args:
            parameter:

        """
        if parameter.number_of_maximum_mismatches == 1:
            df = self.all_mutation(df_input_sequence, parameter)
            # Output as an excel file
            writer = pd.ExcelWriter(f"Mismatched target ({parameter.number_of_maximum_mismatches} bp).xlsx",
                                    engine="xlsxwriter")
            df.to_excel(writer, sheet_name=f"{parameter.number_of_maximum_mismatches} bp MM")
            writer.save()
        else:
            # debug
            # df = self.transversion(df_input_sequence, parameter, 4)

            # 1 mismatches
            df = self.all_mutation(df_input_sequence, parameter)
            writer = pd.ExcelWriter(f"Mismatched target (1 bp).xlsx",
                                    engine="xlsxwriter")
            df.to_excel(writer, sheet_name=f"1 bp MM")
            writer.save()
            # ~ n mismatches
            for i in range(2, parameter.number_of_maximum_mismatches + 1):  # [2,num]
                df = self.transversion(df_input_sequence, parameter, i)
                # Output as an excel file
                writer = pd.ExcelWriter(f"Mismatched target ({i} bp).xlsx",
                                        engine="xlsxwriter")
                df.to_excel(writer, sheet_name=f"{i} bp MM")
                writer.save()


# PAM variants
def PAM(df_input_sequence, parameter: guide_parameter):
    """Generate all possible oligonucleotide sequence with all PAM sequence for length "n" PAM sequence

    making excel spreadsheet

    Input sequence length (format) should be considered carefully.


    Args:
        parameter:


    parameter.
    """
    # Generating all possible PAM sequence
    # First, getting cartesian product of N (PAM length) lists
    pam_string_collection = list()
    pam_obj = [nt for i in range(parameter.len_pam)]
    for combination in itertools.product(*pam_obj):
        string_temp = str()
        for char in combination:
            string_temp += char
        pam_string_collection.append(string_temp)

    # dataframe for PAM data handling
    df_PAM = pd.DataFrame(
        columns=[
            "PAM",
            "Seed guide sequence",
            "PAM + guide",
            "Generated sequence with gDNA context",
            "Gene Name",
            "Known indel frequency from the source"
        ]
    )
    df_PAM.astype({"Known indel frequency from the source": "int32"})
    skipped = int()

    for index, row in df_input_sequence.iterrows():
        for p in pam_string_collection:
            try:
                # TODO: move to input processor
                # input sequence length qualification
                # if len(seed[gene_name]) < SEQUENCE_LENGTH_FOR_VALIDATION:
                #     print("Sequence information is invalid. The operation is aborted.")
                #     return -1

                # input seed pre-processing
                # skip the row with blank cell(s)
                if pd.isna(row["#Sequence with context (150 + 4 + 20 + 126)"]):
                    skipped += 1
                    continue
                leading_sequence, original_pam, protospacer, trailing_sequence = sequence_partitioner(row, parameter)

                df_PAM = df_PAM.append(
                    {
                        "PAM": f"{p}",
                        "Seed guide sequence": f"{protospacer}",
                        "PAM + guide": f"{p}{protospacer}",
                        "Generated sequence with gDNA context": f"{leading_sequence}{p}{protospacer}{trailing_sequence}",
                        "Gene Name": f"{row['#Gene Name']}",
                        "Known indel frequency from the source": f"{row['#Indel frequency']}"
                    },
                    ignore_index=True,
                )
            except Exception as e:
                print(e)
                print("PAM library generation aborted")
                raise

    print(f"{skipped} occurrences are skipped")
    # Output as an excel file
    writer = pd.ExcelWriter("PAM Library.xlsx", engine="xlsxwriter")
    df_PAM.to_excel(writer, sheet_name=f"{len(df_input_sequence.index)} guides processed")
    writer.save()


def guide_random_generator(n, parameter):
    """_summary_

    Args:
        parameter:
    """
    guide_collection = pd.DataFrame(
        columns=[
            "20 bp Randomly generated sequences",
        ]
    )

    # Randomly generating unique guide sequences
    for i in range(int(n)):
        for p in PAM_FOR_RANDOM_GENERTION:
            string_temp = p

            for j in range(parameter.len_guide_seq - 4):
                string_temp += nt[random.sample(range(0, len(nt)), 1)[0]]

            # TODO: is there any better way
            # no duplicates
            if (
                    guide_collection["20 bp Randomly generated sequences"]
                            .str.contains(string_temp)
                            .any()
            ):
                continue
            guide_collection = guide_collection.append(
                {"20 bp Randomly generated sequences": f"{string_temp}"},
                ignore_index=True,
            )
            print(f"{i}-th random sequence has been generated.")

    # Output as an excel file
    writer = pd.ExcelWriter(
        "20 bp Randomly generated sequences.xlsx", engine="xlsxwriter"
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

# Setting the parameter
param = guide_parameter()

src = "input_template.xlsx"
# hard-coded operations
mt = Mutation()
PAM(input_processor(src), param)
mt.mt(input_processor(src), param)
guide_random_generator(NUMBER_OF_FULLY_MATCHED_TARGET, param)
