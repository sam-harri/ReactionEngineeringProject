from Helper.Polynomial import Polynomial
from ReactorEquations import ReactorEquations
from NumericalMethods import NumericalMethods
from ReactorConstants import ReactorConstants
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

#main method
if(__name__ == "__main__"):
    
    # tempArr = np.linspace(600,800)
    # pressureArr = np.linspace(1,10)
    # feedrateArr = np.linspace(1,10)
    # phenolfractionArr = np.linspace(0,0.1)
    
    # for T in temp:
    #     for P in pressure:
    #         for F in feedrate:
    #             for PF in phenolfraction:
    #                 csv.write(method(T,P,F,PF)) #method would return [Volume, Selectivity, ...]
    
    
    #After, normalize each parameter (value/max for things we want to maximize, min/value for things we want to minimize) such that they are bounded [0,1] and 1 being the best
    #Add all normlized parameters per reactor, reactor setting with highest val is the best overall pick
    #Parameters can also be weighted differently with coefficients in front of normalized value (ex : 2*NormV + 1*NormSelec, places 2x more importance on reactor size)
    
    pass