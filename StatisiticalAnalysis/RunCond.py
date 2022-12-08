import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
df = pd.read_csv("CsvData/fullReactorData.csv")

####################################################################################################################################################################
#Data Prep
####################################################################################################################################################################
dfTemp = df[(df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
temp = np.array(dfTemp["Temperature"].values)
tempEff = np.array(dfTemp["Overall Efficiency"].values)

dfPressure = df[(df["Temperature"]==973.15) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
pressure = np.array(dfPressure["Pressure"].values)
pressureEff = np.array(dfPressure["Overall Efficiency"].values)

dfFeed = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
feed = np.array(dfFeed["Feed Rate"].values)
feedEff = np.array(dfFeed["Overall Efficiency"].values)

dfPhenol = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Cross Sectional Area"]==0.2)]
phenol = np.array(dfPhenol["Phenol Fraction"].values)
phenolEff = np.array(dfPhenol["Overall Efficiency"].values)

dfArea = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325)]
area = np.array(dfArea["Cross Sectional Area"].values)
areaEff = np.array(dfArea["Overall Efficiency"].values)

####################################################################################################################################################################
#Model Fitting
####################################################################################################################################################################
tempModel = LinearRegression()
tempModel.fit(temp.reshape(-1,1), tempEff.reshape(-1,1))
tempX = np.linspace(temp.min(), temp.max(), 5)
tempY = tempModel.predict(tempX.reshape(-1,1))

pressureModel = LinearRegression()
pressureModel.fit(pressure.reshape(-1,1), pressureEff.reshape(-1,1))
pressureX = np.linspace(pressure.min(), pressure.max(), 5)
pressureY = pressureModel.predict(pressureX.reshape(-1,1))

feedModel = LinearRegression()
feedModel.fit(feed.reshape(-1,1), feedEff.reshape(-1,1))
feedX = np.linspace(feed.min(), feed.max(), 5)
feedY = feedModel.predict(feedX.reshape(-1,1))

phenolModel = LinearRegression()
phenolModel.fit(phenol.reshape(-1,1), phenolEff.reshape(-1,1))
phenolX = np.linspace(phenol.min(), phenol.max(), 5)
phenolY = phenolModel.predict(phenolX.reshape(-1,1))

areaModel = LinearRegression()
areaModel.fit(area.reshape(-1,1), areaEff.reshape(-1,1))
areaX = np.linspace(area.min(), area.max(), 5)
areaY = areaModel.predict(areaX.reshape(-1,1))

####################################################################################################################################################################
#Plotting Data
####################################################################################################################################################################
plt.subplots_adjust(hspace=0.5, wspace=5)

ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

ax1.scatter(temp, tempEff, c='dodgerblue')
ax1.plot(tempX, tempY, c='dodgerblue')
ax1.set_title("Température Variable")
ax1.set_xlabel("Température [K]")
ax1.set_ylabel("Efficacité")

ax2.scatter(pressure, pressureEff, c='darkcyan')
ax2.plot(pressureX, pressureY, c='darkcyan')
ax2.set_title("Pression Variable")
ax2.set_xlabel("Pression [atm]")
ax2.set_ylabel("Efficacité")

ax3.scatter(feed, feedEff, c='magenta')
ax3.plot(feedX, feedY, c='magenta')
ax3.set_title("Débit Variable")
ax3.set_xlabel("Débit Molaire [mol/s]")
ax3.set_ylabel("Efficacité")

ax4.scatter(phenol, phenolEff, c='firebrick')
ax4.plot(phenolX, phenolY, c='firebrick')
ax4.set_title("Fraction de Phenol Variable")
ax4.set_xlabel("Fraction de Phenol [mol/mol]")
ax4.set_ylabel("Efficacité")

ax5.scatter(area, areaEff, c='darkmagenta')
ax5.plot(areaX, areaY, c='darkmagenta')
ax5.set_ylim(3.121, 3.123)
ax5.set_title("Aire Variable")
ax5.set_xlabel("Aire [m^2]")
ax5.set_ylabel("Efficacité")

plt.show()