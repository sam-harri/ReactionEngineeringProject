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
        Loi de vitesse : water-gas shift (réversible)
        Température en K, concentrations dm^3/L
        """       
        dp : float = -alpha*temp*FT*dW/(2*p*T0*FT0)

        return dp
    
    @staticmethod
    def concentration(P0:float, P:float, T0:float, T:float, F_i:float, F_T:float):
        """
        Calcule la concentration d'une espèce pour une itération
        """
        C_T0 : float = P0/(8.314*T0)

        concentration : float = C_T0*(F_i/F_T)*(P/P0)*(T0/T)

        return concentration

    @staticmethod
    def conversion(F_Ph0:float, F_Ph:float):
        """
        Calcule la conversion du PhOH
        """
        X_Ph : float = (F_Ph0 - F_Ph)/F_Ph0

        return X_Ph

    @staticmethod
    def select_inst(r_H2:float,r_CO:float):
        """
        Calcule la sélectivité instantée du H2 par rapport au CO
        On ne peut le faire par rapport au CO2, car il est produit par la même rxn que l'H2
        """
        S_H2CO : float = r_H2/r_CO

        return S_H2CO

    @staticmethod
    def select_glo(F_H2:float,F_CO:float):
        """
        Calcule la sélectivité globale du H2 par rapport au CO
        On ne peut le faire par rapport au CO2, car il est produit par la même rxn que l'H2
        """
        S_H2CO : float = F_H2/F_CO

        return S_H2CO