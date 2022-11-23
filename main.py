from Helper.Polynomial import Polynomial
from ReactorEquations import ReactorEquations
from NumericalMethods import NumericalMethods
from ReactorConstants import ReactorConstants
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import csv


# #Lois de vitesses
# r_Ph:float = -ReactorEquations.r_srp(temp,C_P,C_H2O)
# r_H2O:float = -5*ReactorEquations.r_srp(temp,C_P,C_H2O)-ReactorEquations.r_wgs(temp,C_P,C_CO,C_H2O,C_CO2,C_H2)
# r_H2:float = 8*ReactorEquations.r_srp(temp,C_P,C_H2O)+ReactorEquations.r_wgs(temp,C_P,C_CO,C_H2O,C_CO2,C_H2)
# r_CO:float = 6*ReactorEquations.r_srp(temp,C_P,C_H2O)-ReactorEquations.r_wgs(temp,C_P,C_CO,C_H2O,C_CO2,C_H2)
# r_CO2:float = ReactorEquations.r_wgs(temp,C_P,C_CO,C_H2O,C_CO2,C_H2)

# #Débit total
# F_T = F_Ph + F_H2O + F_H2 + F_CO + F_CO2

# Calculer "a" (surface extérieure du réacteur/kg catalyseur) selon le A_c (cross-sectionnal area)





#main method
if(__name__ == "__main__"):
    
    csvData = open("rawReactorData.csv", 'w')
    csvData.truncate()
    csvData.close()
    
    csvData = open("processedReactorData.csv", 'w')
    csvData.truncate()
    csvData.close()
    
    rawDF = pd.DataFrame(
        {"Run" : [],
        "Temperature" : [], 
        "Pressure" : [],
        "Feed Rate" : [],
        "Phenol Fraction" : [],
        "Cross Sectional Area" : [],
        "Volume" : [],
        "Selectivity" : [],
        "Conversion" : []}
    )
    
    temperatureRange = np.linspace(600,800, num=5)
    inletPressureRange = np.linspace(1,10, num=5)
    feedRateRange = np.linspace(32,500,num=5)
    phenolFractionRange = np.linspace(0.01, 0.1, num=5)
    crossSectionAreaRange = np.linspace(0.1, 1, num=5)
    
    # testing the addition of rows to the data
    rawDF.loc[len(rawDF.index)] = [0,600,1,30,0.75,0.10,110,12,0.85]
    rawDF.loc[len(rawDF.index)] = [1,650,2,33,1.00,0.05,100,11,0.90]
    
    # run = 0
    # for inletTemperature in temperatureRange:
    #     for inletPressure in inletPressureRange:
    #         for feedRate in feedRateRange:
    #             for phenolFraction in phenolFractionRange:
    #                 for crossSectionalArea in crossSectionAreaRange:
    #                     volume, selectivity, conversion = ReactorEquations.dimentionlizeReactor(inletTemperature,inletPressure,feedRate,phenolFraction, crossSectionalArea)
    #                     rawDF.loc[len(rawDF.index)] = [inletTemperature,
    #                                             inletPressure,
    #                                             feedRate,
    #                                             phenolFraction,
    #                                             crossSectionalArea,
    #                                             volume,
    #                                             selectivity,
    #                                             conversion]
    #                     run +=1
    
    
    rawDF.to_csv("rawReactorData.csv", encoding='utf-8', index=False)
    
    volMin = rawDF["Volume"].min()
    selMax = rawDF["Selectivity"].max()
    convMax = rawDF["Conversion"].max()
    
    runNums = range(len(rawDF.index))
    normVolArr = [volMin/i for i in rawDF["Volume"]]
    normSelArr = [i/selMax for i in rawDF["Selectivity"]]
    normConvArr = [i/convMax for i in rawDF["Conversion"]]
    overallRating = [normVolArr[i]+normSelArr[i]+normConvArr[i] for i in runNums]
    
    normDF = pd.DataFrame({
        "Run" : runNums,
        "Normalized Volume" : normVolArr,
        "Normalized Selectivity" : normSelArr,
        "Normalized Conversion" : normConvArr,
        "Overall Efficiency" : overallRating 
    })
    
    normDF.to_csv("processedReactorData.csv", encoding='utf-8', index=False)
    
    print("Done :)")
    
    
    #After, normalize each parameter (value/max for things we want to maximize, min/value for things we want to minimize) such that they are bounded [0,1] and 1 being the best
    #Add all normlized parameters per reactor, reactor setting with highest val is the best overall pick
    #Parameters can also be weighted differently with coefficients in front of normalized value (ex : 2*NormV + 1*NormSelec, places 2x more importance on reactor size)
