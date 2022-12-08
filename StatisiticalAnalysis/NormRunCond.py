import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
df = pd.read_csv("CsvData/fullReactorData.csv")

####################################################################################################################################################################
#Data Prep
####################################################################################################################################################################
#temp graph Data
dfTemp = df[(df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
temp = np.array(dfTemp["Temperature"].values)
tempEff = np.array(dfTemp["Overall Efficiency"].values)

#pressure graph Data
dfPressure = df[(df["Temperature"]==973.15) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
pressure = np.array(dfPressure["Pressure"].values)
pressureEff = np.array(dfPressure["Overall Efficiency"].values)

#Feed Rate graph
dfFeed = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Phenol Fraction"]==0.0325) & (df["Cross Sectional Area"]==0.2)]
feed = np.array(dfFeed["Feed Rate"].values)
feedEff = np.array(dfFeed["Overall Efficiency"].values)

#Phenol Fraction graph
dfPhenol = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Cross Sectional Area"]==0.2)]
phenol = np.array(dfPhenol["Phenol Fraction"].values)
phenolEff = np.array(dfPhenol["Overall Efficiency"].values)

#Cross Sectional Area graph
dfArea = df[(df["Temperature"]==973.15) & (df["Pressure"]==15.0) & (df["Feed Rate"]==150.0) & (df["Phenol Fraction"]==0.0325)]
area = np.array(dfArea["Cross Sectional Area"].values)
areaEff = np.array(dfArea["Overall Efficiency"].values)

####################################################################################################################################################################
#Normalizing Data
####################################################################################################################################################################
tempMin, tempMax = temp.min(), temp.max()
for i in range(len(temp)):
    temp[i] = (temp[i]-tempMin)/ (tempMax-tempMin)
    
pressureMin, pressureMax = pressure.min(), pressure.max()
for i in range(len(pressure)):
    pressure[i] = (pressure[i]-pressureMin)/ (pressureMax-pressureMin)
    
feedMin, feedMax = feed.min(), feed.max()
for i in range(len(feed)):
    feed[i] = (feed[i]-feedMin)/ (feedMax-feedMin)
    
phenolMin, phenolMax = phenol.min(), phenol.max()
for i in range(len(phenol)):
    phenol[i] = (phenol[i]-phenolMin)/ (phenolMax-phenolMin)
    
areaMin, areaMax = area.min(), area.max()
for i in range(len(area)):
    area[i] = (area[i]-areaMin)/ (areaMax-areaMin)


####################################################################################################################################################################
#Model Building
####################################################################################################################################################################
tempModel = LinearRegression()
tempModel.fit(temp.reshape(-1,1), tempEff.reshape(-1,1))
tempX = np.linspace(0, 1, 5)
tempY = tempModel.predict(tempX.reshape(-1,1))
tempEq = "y = " + str(tempModel.coef_[0][0])[:6] + "x + " + str(tempModel.intercept_[0])[:6]

pressureModel = LinearRegression()
pressureModel.fit(pressure.reshape(-1,1), pressureEff.reshape(-1,1))
pressureX = np.linspace(0, 1, 5)
pressureY = pressureModel.predict(pressureX.reshape(-1,1))
pressureEq = "y = " + str(pressureModel.coef_[0][0])[:6] + "x + " + str(pressureModel.intercept_[0])[:6]

feedModel = LinearRegression()
feedModel.fit(feed.reshape(-1,1), feedEff.reshape(-1,1))
feedX = np.linspace(0, 1, 5)
feedY = feedModel.predict(feedX.reshape(-1,1))
feedEq = "y = " + str(feedModel.coef_[0][0])[:6] + "x + " + str(feedModel.intercept_[0])[:6]

phenolModel = LinearRegression()
phenolModel.fit(phenol.reshape(-1,1), phenolEff.reshape(-1,1))
phenolX = np.linspace(0, 1, 5)
phenolY = phenolModel.predict(phenolX.reshape(-1,1))
phenolEq = "y = " + str(phenolModel.coef_[0][0])[:6] + "x + " + str(phenolModel.intercept_[0])[:6]

areaModel = LinearRegression()
areaModel.fit(area.reshape(-1,1), areaEff.reshape(-1,1))
areaX = np.linspace(0, 1, 5)
areaY = areaModel.predict(areaX.reshape(-1,1))
areaEq = "y = " + str(0.00005619) + "x + " + str(areaModel.intercept_[0])[:6]
#areaModel.coef_[0][0] = 0.00005619

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
ax1.text(.01, .99, tempEq, ha='left', va='top', transform=ax1.transAxes)

ax2.scatter(pressure, pressureEff, c='darkcyan')
ax2.plot(pressureX, pressureY, c='darkcyan')
ax2.set_title("Pression Variable")
ax2.set_xlabel("Pression [atm]")
ax2.set_ylabel("Efficacité")
ax2.text(.01, .99, pressureEq, ha='left', va='top', transform=ax2.transAxes)

ax3.scatter(feed, feedEff, c='magenta')
ax3.plot(feedX, feedY, c='magenta')
ax3.set_title("Débit Variable")
ax3.set_xlabel("Débit Molaire [mol/s]")
ax3.set_ylabel("Efficacité")
ax3.text(.01, .99, feedEq, ha='left', va='top', transform=ax3.transAxes)

ax4.scatter(phenol, phenolEff, c='firebrick')
ax4.plot(phenolX, phenolY, c='firebrick')
ax4.set_title("Fraction de Phenol Variable")
ax4.set_xlabel("Fraction de Phenol [mol/mol]")
ax4.set_ylabel("Efficacité")
ax4.text(.01, .99, phenolEq, ha='left', va='top', transform=ax4.transAxes)

ax5.scatter(area, areaEff, c='darkmagenta')
ax5.plot(areaX, areaY, c='darkmagenta')
ax5.set_ylim(3.121, 3.123)
ax5.set_title("Aire Variable")
ax5.set_xlabel("Aire [m^2]")
ax5.set_ylabel("Efficacité")
ax5.text(.01, .99, areaEq, ha='left', va='top', transform=ax5.transAxes)

plt.show()