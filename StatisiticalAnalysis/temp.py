import pandas as pd

df = pd.read_csv("CsvData/fullReactorData.csv")

df = df.sort_values(by="Overall Efficiency", ascending=False)

cols = ["Weigth","Selectivity","Conversion","Yield","Normalized Weigth","Normalized Selectivity","Normalized Conversion","Normalized Yield"]

for name in cols:
    df.pop(name)

df.to_csv("CsvData/RunCond&Eff.csv", encoding="utf-8", index=False)