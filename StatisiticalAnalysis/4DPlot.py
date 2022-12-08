from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("CsvData2/rawReactorData.csv")

x = np.array(df["Feed Rate"].values)
y = np.array(df["Phenol Fraction"].values)
c = np.array(df["Temperature"].values)
z = np.array(df["Conversion"].values)

cols = ["Run","Temperature","Pressure","Cross Sectional Area","Weigth","Selectivity","Yield"]
for name in cols:
    df.pop(name)
    
df.drop_duplicates(subset=["Feed Rate", "Phenol Fraction"])

feedRate, phenolFraction= np.meshgrid(np.linspace(150,300,num=25), 
                                    np.linspace(0.01, 0.1, num=25))
surfaceX = pd.DataFrame({'FeedRate' : feedRate.ravel(), 'Phenol Fraction' : phenolFraction.ravel()})

reactorModel = make_pipeline(
    PolynomialFeatures(degree=3),
    LinearRegression()
)

reactorModel.fit(df[["Feed Rate", "Phenol Fraction"]],
    df[['Conversion']])

predictedSurface_Conv = reactorModel.predict(surfaceX)

figure = plt.figure(figsize=(15,10))

axis = figure.add_subplot(111, projection='3d')
img = axis.scatter(x, y, z, c=c, cmap=plt.hot())
figure.colorbar(img, label="Temperature [K]")
axis.plot_wireframe(feedRate, phenolFraction, predictedSurface_Conv.reshape(feedRate.shape), alpha=0.3, color='red')
axis.set_title("Reactor Efficiency with Respect to Feed Rate and Phenol Fraction")
axis.set_xlabel("Feed Rate [mol/s]")
axis.set_ylabel("Phenol Fraction [mol/mol]")
axis.set_zlabel("Efficiency")

plt.show()