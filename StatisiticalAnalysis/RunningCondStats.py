import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
df = pd.read_csv("CsvData/fullReactorData.csv")

#Run     Temperature   Pressure   Feed Rate   Phenol Fraction    Cross Sectional Area   ||Weigth,Selectivity,Conversion,Yield
#1759.0  973.15        15.0       150.0       0.0325             0.2                    ||0.0056,1.3333370882082305,0.6442465974217444,8.000009655377054

#temp graph Data
dfTemp = df[(df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
temp = dfTemp["Temperature"]
tempEff = dfTemp["Overall Efficiency"]

#pressure graph Data
dfPressure = df[(df["Temperature"]==973.15) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
pressure = dfPressure["Pressure"]
pressureEff = dfPressure["Overall Efficiency"]

#Feed Rate graph
dfFeed = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
feed = dfFeed["Feed Rate"]
feedArr = dfFeed["Overall Efficiency"]

#Phenol Fraction graph
dfPhenol = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Cross Sectional Area"]==0.2)]
phenol = dfPhenol["Phenol Fraction"]
phenolArr = dfPhenol["Overall Efficiency"]

#Cross Sectional Area graph
dfArea = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325)]
area = dfArea["Cross Sectional Area"]
areaEff = dfArea["Overall Efficiency"]

#Modeling Section
tempModel = LinearRegression()
tempModel.fit(dfTemp[["Temperature"]].values, dfTemp[["Overall Efficiency"]].values)
tempX = np.linspace(dfTemp["Temperature"].min(), dfTemp["Temperature"].max(), 5)
tempY = tempModel.predict(tempX.reshape(-1,1))

pressureModel = LinearRegression()
pressureModel.fit(dfPressure[["Pressure"]].values, dfPressure[["Overall Efficiency"]].values)
pressureX = np.linspace(dfPressure["Pressure"].min(), dfPressure["Pressure"].max(), 5)
pressureY = pressureModel.predict(pressureX.reshape(-1,1))

feedModel = LinearRegression()
feedModel.fit(dfFeed[["Feed Rate"]].values, dfFeed[["Overall Efficiency"]].values)
feedX = np.linspace(dfFeed["Feed Rate"].min(), dfFeed["Feed Rate"].max(), 5)
feedY = feedModel.predict(feedX.reshape(-1,1))

phenolModel = LinearRegression()
phenolModel.fit(dfPhenol[["Phenol Fraction"]].values, dfPhenol[["Overall Efficiency"]].values)
phenolX = np.linspace(dfPhenol["Phenol Fraction"].min(), dfPhenol["Phenol Fraction"].max(), 5)
phenolY = phenolModel.predict(phenolX.reshape(-1,1))

areaModel = LinearRegression()
areaModel.fit(dfArea[["Cross Sectional Area"]].values, dfArea[["Overall Efficiency"]].values)
areaX = np.linspace(dfArea["Cross Sectional Area"].min(), dfArea["Cross Sectional Area"].max(), 5)
areaY = areaModel.predict(areaX.reshape(-1,1))


#Graph Creation
plt.subplots_adjust(hspace=0.5, wspace=5)

ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

ax1.scatter(temp, tempEff, c='dodgerblue')
ax1.plot(tempX, tempY, c='dodgerblue')
ax1.set_title("Variable Temperature")
ax1.set_xlabel("Temperature [K]")
ax1.set_ylabel("Dimentionless Efficiency")

ax2.scatter(pressure, pressureEff, c='crimson')
ax2.plot(pressureX, pressureY, c='crimson')
ax2.set_title("Variable Pressure")
ax2.set_xlabel("Pressure [atm]")
ax2.set_ylabel("Dimentionless Efficiency")

ax3.scatter(feed, feedArr, c='magenta')
ax3.plot(feedX, feedY, c='magenta')
ax3.set_title("Variable Feed")
ax3.set_xlabel("Feed Rate [mol/s]")
ax3.set_ylabel("Dimentionless Efficiency")

ax4.scatter(phenol, phenolArr, c='mediumorchid')
ax4.plot(phenolX, phenolY, c='mediumorchid')
ax4.set_title("Variable Phenol Fraction")
ax4.set_xlabel("Phenol Fraction [mol/mol]")
ax4.set_ylabel("Dimentionless Efficiency")

ax5.scatter(area, areaEff, c='darkorange')
ax5.plot(areaX, areaY, c='darkorange')
ax5.set_ylim(3.121, 3.123)
ax5.set_title("Variable Area")
ax5.set_xlabel("Cross Sectional Area [m^2]")
ax5.set_ylabel("Dimentionless Efficiency")

plt.show()