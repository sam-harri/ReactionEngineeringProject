from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from DimentionReactor import DimentionReactor
import numpy as np
import pandas as pd
import time

start = time.time()

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
        "Weigth" : [],
        "Selectivity" : [],
        "Conversion" : [],
        "Yield" : []}
    )
    
    temperatureRange = np.linspace(873.15,1073.15, num=5)
    inletPressureRange = np.linspace(2,10, num=5)
    feedRateRange = np.linspace(150,300,num=5)
    phenolFractionRange = np.linspace(0.01, 0.1, num=5)
    crossSectionAreaRange = np.linspace(0.05, 0.2, num=5)
    
    run = 0
    for inletTemperature in temperatureRange:
        for inletPressure in inletPressureRange:
            for feedRate in feedRateRange:
                for phenolFraction in phenolFractionRange:
                    for crossSectionalArea in crossSectionAreaRange:
                        weight, conversion, selectivity, rendement = DimentionReactor.dimentionlizeReactor(inletTemperature,inletPressure,feedRate,phenolFraction,crossSectionalArea)
                        print(run)
                        tmp = np.array([
                            run,
                            inletTemperature,
                            inletPressure,
                            feedRate,
                            phenolFraction,
                            crossSectionalArea,
                            weight,
                            selectivity,
                            conversion,
                            rendement])
                        rawDF.loc[len(rawDF.index)] = tmp
                        run +=1

    print(rawDF.shape)                    
    index_names = rawDF[rawDF['Weigth'] == 0].index
    rawDF.drop(index_names, inplace = True)
    print(rawDF.shape)
    
    rawDF.to_csv("CsvData/rawReactorData.csv", encoding='utf-8', index=False)
    
    weightMin = rawDF["Weigth"].min()
    selMax = rawDF["Selectivity"].max()
    convMax = rawDF["Conversion"].max()
    rendementMax = rawDF["Yield"].max()
    
    runNums = range(len(rawDF.index))
    normWeiArr = [weightMin/i for i in rawDF["Weigth"]]
    normSelArr = [i/selMax for i in rawDF["Selectivity"]]
    normConvArr = [i/convMax for i in rawDF["Conversion"]]
    normRendArr = [i/rendementMax for i in rawDF["Yield"]]
    overallRating = [normWeiArr[i]+normSelArr[i]+normConvArr[i]+normRendArr[i] for i in runNums]
    
    normDF = pd.DataFrame({
        "Run" : runNums,
        "Normalized Weigth" : normWeiArr,
        "Normalized Selectivity" : normSelArr,
        "Normalized Conversion" : normConvArr,
        "Normalized Yield" : normRendArr,
        "Overall Efficiency" : overallRating 
    })
    normDF.to_csv("CsvData/processedReactorData.csv", encoding='utf-8', index=False)
    
    fullDF = pd.merge(rawDF, normDF, on="Run", how="inner")
    fullDF.to_csv("CsvData/fullReactorData.csv", encoding="utf-8", index=False)
    
    stop = time.time()
    
    print(str((stop-start)//60) + "minutes and " + str((stop-start)%60) + " seconds")
    
    # reactorModel = make_pipeline(
    #     PolynomialFeatures(degree=3),
    #     LinearRegression()
    # )
    
    # reactorModel.fit(fullDF[[
    #     "Temperature",
    #     "Pressure",
    #     "Feed Rate",
    #     "Phenol Fraction",
    #     "Cross Sectional Area"]],
    #     fullDF[['Overall Efficiency']])
    
    
    
    # print(fullDF)
    
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
    #https://towardsdatascience.com/polynomial-regression-gradient-descent-from-scratch-279db2936fe9
