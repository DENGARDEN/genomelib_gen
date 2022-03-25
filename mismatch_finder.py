import itertools
import random
import re

import pandas as pd


NONWOBBLE_TRANSITION = {"G": ["a"], "C": ["t"], "T": ["c"]}
WOBBLE_TRANSITION = {"A": ["g"]}
TRANSVERSION = {"A": ["c", "t"], "G": ["c", "t"], "C": ["a", "g"], "T": ["a", "g"]}
# TRANSVERSION_SIM_GUIDE = {"A": "t", "G": "c", "C": "g", "T": "a"}


def mismatch_finder(original_sequence: str, mutant_sequence: str):

    mt_pos = str()
    mt_cnt = int()
    mt_type = str()
    
    for idx in range(len(mutant_sequence)):
        if mutant_sequence[idx].islower():
            mt_pos+=f"{idx}, "
            mt_cnt += 1

            if mutant_sequence[idx] in TRANSVERSION[original_sequence[idx]]  :
                mt_type+=f"TRANSVERSION, "
            elif original_sequence[idx] == "A":
                if mutant_sequence[idx] in WOBBLE_TRANSITION[original_sequence[idx]]  :
                    mt_type+=f"WOBBLE_TRANSITION, "
            elif mutant_sequence[idx] in NONWOBBLE_TRANSITION[original_sequence[idx]]  :
                mt_type+=f"NONWOBBLE_TRANSITION, "

    return mt_pos, mt_cnt, mt_type


src = "2017_NatMethods_Kim_highthroughput Cpf1_mismatch guide.xlsx"
df = pd.read_excel(src, engine="openpyxl", skiprows=[0])
print(df)
df["Mismatch position"]= df["Mismatch position"].astype("string")
df["Mutation type"]= df["Mutation type"].astype("string")

for index, row in df.iterrows():
    pos, n, type = mismatch_finder(
        row["Guide sequence (5' to 3')"], row["Target sequence (5' to 3')"]
    )
    df.at[index, "Mismatch position"] = pos
    df.at[index, "Number of mismatched bases (bp)"] = n
    df.at[index, "Mutation type"] = type

writer = pd.ExcelWriter(f"hidden data revealed.xlsx", engine="xlsxwriter")
df.to_excel(writer, sheet_name=f"attributes added")
writer.save()
