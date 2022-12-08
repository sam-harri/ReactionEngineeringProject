import pandas as pd

df = pd.read_csv("CsvData2/rawReactorData.csv")
df = df.sort_values(by="Conversion", ascending=False)
df.to_csv("CsvData2/bestConversion.csv")
df = df.sort_values(by="Weigth", ascending=True)
df.to_csv("CsvData2/bestWeigth.csv")