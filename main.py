# from Helper.Polynomial import Polynomial
# from ReactorEquations import ReactorEquations
# from NumericalMethods import NumericalMethods
# from ReactorConstants import ReactorConstants
# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd
# import math


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
    pass

    # tempArr = np.linspace(600,800)
    # pressureArr = np.linspace(1,10)
    # feedrateArr = np.linspace(1,10)
    # phenolfractionArr = np.linspace(0,0.1)
    


#USE A WHILE LOOP UNTIL WE CONVERGE TO F_H2 = 25 MOL/S ??
#WE SHOULD HAVE A FOR LOOP FOR THE CROSS SECTIONAL AREA



    # for T in tempArr:
    #    for P in pressureArr:
    #       for F in feedrateArr:
    #            for PF in phenolfractionArr:              
    #                 csv.write(method(T,P,F,PF)) #method would return [Volume, Selectivity,...]
    
    
    #After, normalize each parameter (value/max for things we want to maximize, min/value for things we want to minimize) such that they are bounded [0,1] and 1 being the best
    #Add all normlized parameters per reactor, reactor setting with highest val is the best overall pick
    #Parameters can also be weighted differently with coefficients in front of normalized value (ex : 2*NormV + 1*NormSelec, places 2x more importance on reactor size)
