import matplotlib.pyplot as plt
import numpy as np
from Helper.Polynomial import Polynomial
from ReactorConstants import ReactorConstants
from typing import List
from math import sqrt,exp

F_H2 : float = ReactorConstants.F_H2
K_P : Polynomial = ReactorConstants.KP
K_H2O : Polynomial = ReactorConstants.KH2O
D_eq : float = ReactorConstants.D_eq
phi : float = ReactorConstants.phi
rho_c : float = ReactorConstants.rho_c
U : float = ReactorConstants.U
mu: float = ReactorConstants.mu
Cp_Ph_model : Polynomial = ReactorConstants.Cp_Ph_model
Cp_calo : float = ReactorConstants.Cp_calo

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
    def alpha(P:float,temp:float, x_Ph:float, A_c:float,P0:float,FT:float):
        """
        Alpha de l'équation d'ergun
        T en K et P en Pa
        # """
        # rho_0 est la masse volumique initiale du gaz
        # rho_Ph est la masse volumique du PhOH au conditions initiales
        # rho_c est la masse volumique des particules de catalyseur
        débit_vol : float = FT*8.314*temp/P
        rho_Ph : float = 0.09411*P/(8.314*temp) # 0.09411 kg/mol, masse molaire du PhOH
        w_Ph : float = x_Ph*94.11/(x_Ph*94.11+(0.9-x_Ph)*18.02+1.401) # x_Ph est la fraction molaire de Ph
        rho_0 : float = rho_Ph/w_Ph
        u : float = débit_vol/A_c
        G : float = rho_0*u

        beta : float = -G*(1-phi)/(rho_0*1*D_eq*phi**3)*(150*(1-phi)*mu/D_eq+1.75*G)
        alpha : float = 2*beta/((1-phi)*A_c*rho_c*P0)

        return alpha
    
    @staticmethod
    def ergun(alpha:float, p:float, temp:float, T_0:float, F_T:float, F_T0:float):
        """
        Loi de vitesse : water-gas shift (réversible)
        Température en K, concentrations dm^3/L
        """       
        dpdW : float = -alpha*temp*F_T/(2*p*T_0*F_T0)

        return dpdW
    
    @staticmethod
    def conception(r_i:float):
        """
        Équation de conception pour un PBR
        """
        return r_i

    @staticmethod
    def bilan_E(T:float,Ta:float,a:float,F_Ph:float,F_H2O:float,F_H2:float,F_CO:float,F_CO2:float,F_N2:float,C_P:float,C_H2O:float,C_CO:float,C_CO2:float,C_H2:float):
        """
        Bilan énergétique pour un PBR à réactions multiples (Équation T11-1.J du livre Folger (2018))
        T en K
        Les corrélations et les coefficients pour les Cp viennent du livre Felder, Rousseau, Bullard (2019)
        Le "a" est la surface d'échange de chaleur par kg de catalyseur
        """
        T_R = 298.15
        K = 273.15
        
        # Enthalpies de formation à 25 deg C et 1 atm en phase gazeuse (Source: Felder et Rousseau)
        H_Ph : float = -90800 + Polynomial.definiteIntegral(Cp_Ph_model,T_R, T) # J/mol
        H_H2O: float = -241830 + 1000*(33.46*10**(-3)*(T-T_R)+0.6880*10**(-5)*1/2*(T-T_R)**2+0.7604*10**(-8)*1/3*(T-T_R)**3-3.593*10**(-12)*1/4*(T-T_R)**4) # J/mol
        H_H2 : float = 0 + 1000*(28.84*10**(-3)*(T-T_R)+0.00765*10**(-5)*1/2*(T-T_R)**2+0.3288*10**(-8)*1/3*(T-T_R)**3-0.8698*10**(-12)*1/4*(T-T_R)**4) # J/mol
        H_CO : float = -110520 + 1000*(28.95*10**(-3)*(T-T_R)+0.4110*10**(-5)*1/2*(T-T_R)**2+0.3548*10**(-8)*1/3*(T-T_R)**3-2.220*10**(-12)*1/4*(T-T_R)**4) # J/mol
        H_CO2: float = -393500 + 1000*(36.11*10**(-3)*(T-T_R)+4.233*10**(-5)*1/2*(T-T_R)**2-2.887*10**(-8)*1/3*(T-T_R)**3+7.464*10**(-12)*1/4*(T-T_R)**4) # J/mol

        deltaH_1 : float = 6/1*H_CO + 8/1*H_H2 - 5/1*H_H2O - H_Ph
        deltaH_2 : float = 1/1*H_CO2 + 1/1*H_H2 - 1/1*H_H2O - H_CO
        
        # Attention, certaines formules ont des unités de K et d'autres de deg C
        Cp_Ph : Polynomial = Cp_Ph_model.evaluate(T) # https://webbook.nist.gov/cgi/cbook.cgi?ID=C108952&Mask=1E9F#Thermo-Gas
        Cp_H2O: float = 1000*(33.46*10**(-3)+0.6880*10**(-5)*(T+K)+0.7604*10**(-8)*(T+K)**2-3.593*10**(-12)*(T+K)**3) #J/(mol*deg C)
        Cp_H2 : float = 1000*(28.84*10**(-3)+0.00765*10**(-5)*(T+K)+0.3288*10**(-8)*(T+K)**2-0.8698*10**(-12)*(T+K)**3) #J/(mol*deg C)
        Cp_CO : float = 1000*(28.95*10**(-3)+0.4110*10**(-5)*(T+K)+0.3548*10**(-8)*(T+K)**2-2.220*10**(-12)*(T+K)**3) #J/(mol*deg C)
        Cp_CO2: float = 1000*(36.11*10**(-3)+4.233*10**(-5)*(T+K)-2.887*10**(-8)*(T+K)**2+7.464*10**(-12)*(T+K)**3) #J/(mol*deg C)
        Cp_N2 : float = 1000*(29.00*10**(-3)+0.2199*10**(-5)*(T+K)+0.5723*10**(-8)*(T+K)**2-2.871*10**(-12)*(T+K)**3) #J/(mol*deg C)

        somme_rxn : float = ReactorEquations.r_srp(T,C_P,C_H2O)*deltaH_1 + ReactorEquations.r_wgs(T,C_P,C_CO,C_H2O,C_CO2,C_H2)*deltaH_2
        somme_debit : float = F_Ph*Cp_Ph + F_H2O*Cp_H2O + F_H2*Cp_H2 + F_CO*Cp_CO + F_CO2*Cp_CO2 + F_N2*Cp_N2

        dTdW : float = (somme_rxn-U*a*(T-Ta))/somme_debit

        return dTdW
    
    @staticmethod
    def échange_chal(T:float,Ta:float,a:float,m:float):
        """
        Équation différentielle pour le fluide caloporteur
        Le "a" est la surface d'échange de chaleur par kg de catalyseur
        Le "m" est le débit massique de fluide caloporteur
        """

        dTadW : float = U*a*(T-Ta)/(m*Cp_calo)

        return dTadW

    @staticmethod
    def concentration(P_0:float, P:float, T_0:float, T:float, F_i:float, F_T:float):
        """
        Calcule la concentration d'une espèce pour une itération
        """
        return (P_0/(8.314*T_0))*(F_i/F_T)*(P/P_0)*(T_0/T)

    @staticmethod
    def conversion(F_Ph0:float, F_Ph:float):
        """
        Calcule la conversion du PhOH
        """
        return (F_Ph0 - F_Ph)/F_Ph0

    @staticmethod
    def select_inst(r_H2:float,r_CO:float):
        """
        Calcule la sélectivité instantée du H2 par rapport au CO
        On ne peut le faire par rapport au CO2, car il est produit par la même rxn que l'H2
        """
        return r_H2/r_CO

    @staticmethod
    def select_glo(F_H2:float,F_CO:float):
        """
        Calcule la sélectivité globale du H2 par rapport au CO
        On ne peut le faire par rapport au CO2, car il est produit par la même rxn que l'H2
        """
        return F_H2/F_CO
    
    @staticmethod
    def dimentionlizeReactor(temperature: float, inletPressure: float, feedRate: float, phenolFraction: float):
        """
        returns volume, selectivity, and conversion
        """
        pass