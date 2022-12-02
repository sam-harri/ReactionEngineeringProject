from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from ReactorEquations import ReactorEquations
import numpy as np
import pandas as pd
import csv

#main method
if(__name__ == "__main__"):
    
    csvFiles = [
        "CsvData/rawReactorData.csv",
        "CsvData/processedReactorData.csv",
        "CsvData/fullReactorData.csv",
        "CsvData/tempRawReactorData.csv"]
    
    for csvData in csvFiles:
        tmp = open(csvData, "w")
        tmp.truncate()
        tmp.close()
    
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
    inletPressureRange = np.linspace(1,5, num=5)
    feedRateRange = np.linspace(32,500,num=5)
    phenolFractionRange = np.linspace(0.01, 0.1, num=5)
    crossSectionAreaRange = np.linspace(0.1, 1, num=5)
    
    
    # with open("CsvData/tempRawReactorData.csv", "wb") as tempData:
    #     tempData.write(",".join([0,600,1,30,0.75,0.10,110,12,0.85]))
    #     tempData.write("\n")
    #     tempData.write(",".join([1,650,2,33,1.00,0.05,100,11,0.90]))
    #     tempData.write("\n")
    # # testing the addition of rows to the data
    # a = [0,600,1,30,0.75,0.10,110,12,0.85]
    # rawDF.loc[len(rawDF.index)] = a
    # b = [1,650,2,33,1.00,0.05,100,11,0.90]
    # rawDF.loc[len(rawDF.index)] = b
    
    run = 0
    for inletTemperature in temperatureRange:
        for inletPressure in inletPressureRange:
            for feedRate in feedRateRange:
                for phenolFraction in phenolFractionRange:
                    for crossSectionalArea in crossSectionAreaRange:
                        volume, selectivity, conversion = ReactorEquations.dimentionlizeReactor(inletTemperature,inletPressure,feedRate,phenolFraction,crossSectionalArea)
                        if(volume==0):
                            print(f'Run {run} : Failure')
                        else:
                            print(f'Run {run} : Success - {volume}, {selectivity}, {conversion}')
                        tmp = np.array([
                            run,
                            inletTemperature,
                            inletPressure,
                            feedRate,
                            phenolFraction,
                            crossSectionalArea,
                            volume,
                            selectivity,
                            conversion])
                        np.savetxt('myfile.csv', tmp, delimiter=',')
                        rawDF.loc[len(rawDF.index)] = tmp
                        run +=1
    
    rawDF.to_csv("CsvData/rawReactorData.csv", encoding='utf-8', index=False)
    
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
    normDF.to_csv("CsvData/processedReactorData.csv", encoding='utf-8', index=False)
    
    fullDF = pd.merge(rawDF, normDF, on="Run", how="inner")
    fullDF.to_csv("CsvData/fullReactorData.csv", encoding="utf-8", index=False)
    
    reactorModel = make_pipeline(
        PolynomialFeatures(degree=3),
        LinearRegression()
    )
    
    reactorModel.fit(fullDF[[
        "Temperature",
        "Pressure",
        "Feed Rate",
        "Phenol Fraction",
        "Cross Sectional Area"]],
        fullDF[['Overall Efficiency']])
    
    
    
    print(fullDF)
    
    # t_surf, p_surf, fr_surf, pf_surf, csa_surf = np.meshgrid(
    #     temperatureRange,
    #     inletPressureRange,
    #     feedRateRange,
    #     phenolFractionRange,
    #     crossSectionAreaRange,
    # )
    
    # surfaceReactor = pd.DataFrame({
    #     "Temperature" : t_surf.ravel(),
    #     "Pressure" : p_surf.ravel(),
    #     "Feed Rate" : fr_surf.ravel(),
    #     "Phenol Fraction" : fr_surf.ravel(),
    #     "Cross Sectional Area" : csa_surf.ravel()
    # })
    
    # predictedReactorData = reactorModel.fit(surfaceReactor)
    
    #https://stackoverflow.com/questions/54859865/scikit-learn-order-of-coefficients-for-multiple-linear-regression-and-polynomia
