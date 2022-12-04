import pandas as pd

df = pd.read_csv("CsvData/fullReactorData.csv")

df = df.sort_values(by="Overall Efficiency", ascending=False)

df.to_csv("CsvData/testFullData.csv", encoding="utf-8", index=False)