import tabula
import pandas as pd

# src = input("type the file name with extension: ")
# page = input("type the page information you want to get: ")

# debug
df = tabula.read_pdf("41589_2021_868_MOESM1_ESM.pdf", pages = "16-17")

# df = tabula.read_pdf(f"{src}", pages = f"{page}")

res = pd.DataFrame()
for i in range(len(df)):
    res =res.append(df[i],ignore_index =True)

res.to_excel("extracted.xlsx")
# res.to_excel(f'{src}_extracted.xlsx')

print(res)
