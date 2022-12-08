import pandas as pd

df = pd.read_csv("CsvData/fullReactorData.csv")

df = df.sort_values(by="Weigth", ascending=True)

df.to_csv("CsvData/bestWeigth.csv", encoding="utf-8", index=False)