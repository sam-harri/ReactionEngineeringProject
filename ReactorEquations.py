import matplotlib.pyplot as plt
import numpy as np
from Helper.Polynomial import Polynomial
from ReactorConstants import ReactorConstants
from typing import List
from math import sqrt,exp

FA0: float = ReactorConstants.FA0
weigth: float = ReactorConstants.weight
K_P : Polynomial = ReactorConstants.KP
K_H2O : Polynomial = ReactorConstants.KH2O
D_eq : float = ReactorConstants.D_eq
phi : float = ReactorConstants.phi
rho_c : float = ReactorConstants.rho_c

class ReactorEquations:

    @staticmethod
    def r_srp(temp: float, C_P: float, C_H2O: float):
        """
        (float,float,float->float)
        Loi de vitesse : reformage du phénol
        Température en K, concentrations dm^3/L
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
        (float,float,float,float,float,float)->float)
        Loi de vitesse : water-gas shift (réversible)
        Température en K, concentrations dm^3/L
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
    def alpha(P:float,temp:float, y_Ph:float, mu:float,A_c:float,P0:float,débit_vol:float):
        """
        (float,float,float,float,float,float,float)->float)
        Alpha de l'équation d'ergun
        rho_0 est la masse volumique initiale du gaz
        rho_Ph est la masse volumique du PhOH au conditions initiales
        rho_c est la masse volumique des particules de catalyseur
        """

        rho_Ph : float = 0.09411*P/(8.314*temp) # 0.09411 kg/mol, masse molaire du PhOH
        x_Ph : float = 94.11/(18.02/y_Ph-18.02+94.11) # y_Ph est la fraction molaire de Ph
        rho_0 : float = rho_Ph/x_Ph
        u : float = débit_vol/A_c
        G : float = rho_0*u

        beta : float = -G*(1-phi)/(rho_0*1*D_eq*phi**3)*(150*(1-phi)*mu/D_eq+1.75*G)
        alpha : float = 2*beta/((1-phi)*A_c*rho_c*P0)

        return alpha
    
    @staticmethod
    def ergun(alpha:float, p:float, temp:float, T0:float, FT:float, FT0:float, dW:float):
        """
        (float,float,float,float,float,float)->float)
        Loi de vitesse : water-gas shift (réversible)
        Température en K, concentrations dm^3/L
        """       
        dp : float = -alpha*temp*FT*dW/(2*p*T0*FT0)

        return dp