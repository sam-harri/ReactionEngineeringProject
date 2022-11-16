import matplotlib.pyplot as plt
import numpy as np
from Helper.Polynomial import Polynomial
from ReactorConstants import ReactorConstants
from typing import List
from math import sqrt,exp

FA0: float = ReactorConstants.FA0
alpha: float = ReactorConstants.alpha
weigth: float = ReactorConstants.weight
K_P : Polynomial = ReactorConstants.KP
K_H2O : Polynomial = ReactorConstants.KH2O

class ReactorEquations:

    @staticmethod
    def r_srp(temp: float, C_P: float, C_H2O: float):
        """
        temp in K
        """
        k_srp: float = exp(7.63-39960/(8.314*temp))
        K_P: float = K_P.evaluate(temp)
        K_H2O: float = K_H2O.evaluate(temp)

        P_P: float = C_P*8.314*temp
        P_H2O: float = C_H2O*8.314*temp

        r_srp: float = k_srp*K_P*K_H2O*(P_P*P_H2O**5)/(K_P*P_P+K_H2O*P_H2O+1)**2

        return r_srp

    @staticmethod
    def r_wgs(temp: float, C_P: float, C_CO: float, C_H2O: float, C_CO2: float, C_H2: float):
        """
        temp in K
        """
        k_wgs: float = exp(7.45-40100/(8.314*temp))
        k_rwgs: float = exp(7.03-38840/(8.314*temp))
        K_wgs: float = k_wgs/k_rwgs
        K_P: float = K_P.evaluate(temp)
        K_H2O: float = K_H2O.evaluate(temp)
        
        P_P: float = C_P*8.314*temp
        P_CO: float = C_CO*8.314*temp
        P_H2O: float = C_H2O*8.314*temp
        P_CO2: float = C_CO2*8.314*temp
        P_H2: float = C_H2*8.314*temp

        r_wgs: float = k_wgs*(P_CO*P_H2O-P_CO2*P_H2/K_wgs)/(K_P*P_P+K_H2O*P_H2O+1)**2

        return r_wgs
    
    @staticmethod
    def ergun(alpha:float, p:float, temp:float, T0:float, FT:float, FT0:float, dW:float):
        """
        p is adimensionnal pressure P/P0
        """       
        return -alpha*T0*FT*dW/(2*p*T0*FT0)