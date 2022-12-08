import matplotlib.pyplot as plt
from typing import List
from math import pi, exp

class Term:
    def __init__(self,coefficient: float, power: int) -> None:
        self.coefficient = coefficient
        self.power = power
        
    def evaluate(self,xval: float):
        return self.coefficient * (pow(xval,self.power))
    
class Polynomial:
    def __init__(self, polyRep: str) -> None:
        self.polyRep : str = polyRep
        self.termArr : List[Term] = Polynomial.toPolynomial(self.polyRep)
        
    @staticmethod
    def toPolynomial(polyRep: str):
        stringArr: List[str] = polyRep.split(";")
        termArr: List[Term] = []
        for string in stringArr:
            coeff, power = string.split("x^")
            termArr.append(Term(float(coeff),int(power)))
        return termArr
    
    def evaluate(self, xval: float):
        sum: float = 0
        for term in self.termArr:
            sum += term.evaluate(xval)
        return sum
    
    def derive(self):
        for term in self.termArr:
            term.coefficient *= term.power
            term.power -= 1
    
    @staticmethod
    def integrate(poly):
        for term in poly.termArr:
            term.coefficient /= term.power+1
            term.power += 1
        return poly
    
    @staticmethod
    def definiteIntegral(poly, lower, upper):
        for term in poly.termArr:
            term.coefficient /= term.power+1
            term.power += 1
        return poly.evaluate(upper) - poly.evaluate(lower)

class ReactorConstants:
    
    #add all constants here so that we can easily update them acrossall classes
    #changes here will be reflected on everything else
    #make sure to delcare these as the constants by using the class name everywhere else though
    #ex : FA0: float = ReactorConstants.FA0
    
    F_H2: float = 25 #molH2/s à la sortie
    YI0: float = 0.1 #molI/mol
    U: float = 16.0 #W/m**2K
    KP: Polynomial = Polynomial("0.000041x^2;-0.0600983x^1;22.88680282250014x^0")
    KH2O: Polynomial = Polynomial("0.0000745x^2;-0.11574935x^1;50.048421226250284x^0")
    Cp_Ph_model : Polynomial = Polynomial("-0.0000782500x^2;0.242785x^1;68.01049999999961x^0")
    CP_Ph_model_integrated : Polynomial = Polynomial.integrate(Cp_Ph_model)
    D_eq: float = 0.005 #m, diamètre équivalent catalyseur
    phi: float = 0.39 #porosité
    rho_c: float = 5740 #kg/m^3, masse volumique catalyseur
    mu: float = 24.59*10**(-6) #Pa*s, viscosité à 700 deg C, 1 atm, 10% mol N2 et 90% mol H2O
    Cp_calo: float = 4187 # J/(kg*K), Cp de l'eau à 15 deg C (https://www.engineeringtoolbox.com/water-thermal-properties)
    T_calo_in: float = 288.15 # K, température de l'eau à l'entrée
    StepSize : float = 0.0001 #kg of catalyst
    
    m_calo: float = 2 #kg/s

F_H2 : float = ReactorConstants.F_H2
K_P : Polynomial = ReactorConstants.KP
K_H2O : Polynomial = ReactorConstants.KH2O
D_eq : float = ReactorConstants.D_eq
phi : float = ReactorConstants.phi
rho_c : float = ReactorConstants.rho_c
U : float = ReactorConstants.U
mu: float = ReactorConstants.mu
Cp_Ph_model : Polynomial = ReactorConstants.Cp_Ph_model
Cp_Ph_model_integrated : Polynomial = ReactorConstants.CP_Ph_model_integrated
Cp_calo : float = ReactorConstants.Cp_calo
YI0: float = ReactorConstants.YI0
Ta0: float = ReactorConstants.T_calo_in
m_calo: float = ReactorConstants.m_calo
h = ReactorConstants.StepSize

class ReactorEquations:

    @staticmethod
    def r_srp(temp: float, C_Ph: float, C_H2O: float):
        """
        temperature in K
        C_j in mol/m^3
        """
        k_srp: float = exp(7.63-39960/(8.314*temp)) / 3.6 # mol / kg * s
        K_P_eval: float = K_P.evaluate(temp) #1 / atm
        K_H2O_eval: float = K_H2O.evaluate(temp) #1 / atm

        P_P: float = C_Ph*8.2057*10**(-5)*temp #atm
        P_H2O: float = C_H2O*8.2057*10**(-5)*temp #atm

        return k_srp*K_P_eval*K_H2O_eval*(P_P*P_H2O**5)/(K_P_eval*P_P+K_H2O_eval*P_H2O+1)**2

    @staticmethod
    def r_wgs(temp: float, C_Ph: float, C_CO: float, C_H2O: float, C_CO2: float, C_H2: float):
        """
        Loi de vitesse : water-gas shift (réversible)
        Température en K, concentrations dm^3/L
        """
        k_wgs: float = exp(7.45-40100/(8.314*temp)) / 3.6 # mol / kg * s
        k_rwgs: float = exp(7.03-38840/(8.314*temp)) / 3.6 # mol / kg * s
        K_wgs: float = k_wgs/k_rwgs #adim
        K_P_eval: float = K_P.evaluate(temp) #1 / atm
        K_H2O_eval: float = K_H2O.evaluate(temp) #1 / atm
        
        P_P: float = C_Ph*8.2057*10**(-5)*temp  # atm
        P_CO: float = C_CO*8.2057*10**(-5)*temp # atm
        P_H2O: float = C_H2O*8.2057*10**(-5)*temp # atm
        P_CO2: float = C_CO2*8.2057*10**(-5)*temp # atm
        P_H2: float = C_H2*8.2057*10**(-5)*temp # atm

        return k_wgs*(P_CO*P_H2O-P_CO2*P_H2/K_wgs)/(K_P_eval*P_P+K_H2O_eval*P_H2O+1)**2

    # Lois de vitesses globales
    @staticmethod
    def r_Ph(r_srp:float):
        return -r_srp

    @staticmethod
    def r_H2O(r_srp:float,r_wgs:float):
        return -5*r_srp - r_wgs

    @staticmethod
    def r_H2(r_srp:float,r_wgs:float):
        return 8*r_srp + r_wgs

    @staticmethod
    def r_CO(r_srp:float,r_wgs:float):
        return 6*r_srp - r_wgs

    @staticmethod
    def r_CO2(r_wgs:float):
        return r_wgs
    
    @staticmethod
    def alpha(P:float,T0:float, temp:float, x_Ph:float, A_c:float,P0:float,FT:float, phenolFraction):
        """
        Alpha de l'équation d'ergun
        T0 en K et P0 en Pa (conditions initiales)
        temp et P sont les valeurs à l'itération
        x_Ph est la fraction molaire initiale de Ph
        """
        # rho_0 est la masse volumique initiale du gaz
        # rho_Ph est la masse volumique du PhOH aux conditions initiales
        # rho_c est la masse volumique des particules de catalyseur
        débit_vol : float = FT*8.2057*10**(-5)*temp/P
        rho_Ph : float = 0.09411*P0/(8.2057*10**(-5)*T0)*phenolFraction # 0.09411 kg/mol, masse molaire du PhOH
        w_Ph : float = x_Ph*94.11/(x_Ph*94.11+(0.9-x_Ph)*18.02+(0.1)*28.014) # w_Ph est la fraction massique initiale de Ph
        rho_0 : float = rho_Ph/w_Ph
        u : float = débit_vol/A_c
        G : float = rho_0*u

        beta : float = (((G*(1-phi))/(rho_0*1*D_eq*phi**3))*(((150*(1-phi)*mu)/D_eq)+(1.75*G))) / 101325
        return 2*beta/((1-phi)*A_c*rho_c*P0)

    @staticmethod
    def ergun(alpha:float, p:float, temp:float, T_0:float, F_T:float, F_T0:float):
        """
        Eq chute de pression
        """       
        return -alpha*temp*F_T/(2*p*T_0*F_T0)
    
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
        K = -273.15
        
        # Enthalpies de formation à 25 deg C et 1 atm en phase gazeuse (Source: Felder et Rousseau)
        H_Ph : float = -90800 + Cp_Ph_model_integrated.evaluate(T)-Cp_Ph_model_integrated.evaluate(T_R) # J/mol
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

        dTdW : float = (-somme_rxn-U*a*(T-Ta))/somme_debit

        return dTdW

    @staticmethod
    def echange_chal(T:float,Ta:float,a:float,m:float):
        """
        Équation différentielle pour le fluide caloporteur
        Le "a" est la surface d'échange de chaleur par kg de catalyseur
        Le "m" est le débit massique de fluide caloporteur
        """

        dTadW : float = U*a*(T-Ta)/(m*Cp_calo)

        return dTadW

    @staticmethod
    def a(Ac:float):
        """
        Surface d'échange de chaleur par kg de catalyseur
        """

        D : float = Ac*4/pi
        rho_bulk : float = rho_c*(1-phi)
        a : float = 4/(D*rho_bulk)

        return a
    
    @staticmethod
    def concentration(P_0:float, P:float, T_0:float, T:float, F_i:float, F_T:float):
        """
        Calcule la concentration d'une espèce pour une itération
        """
        return (P_0/(8.2057*10**(-5)*T_0))*(F_i/F_T)*(P/P_0)*(T_0/T)

    @staticmethod
    def F_T(F_Ph:float,F_H2O:float,F_H_2:float,F_CO:float,F_CO2:float, F_N2: float):
        """
        Débit total
        """
        F_T : float = F_Ph + F_H2O + F_H_2 + F_CO + F_CO2 + F_N2

        return F_T

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
    def rendement(F_Ph0:float,F_Ph:float,F_H2:float):
        """
        Calcule le rendement globale du H2 par rapport au Ph ayant réagit
        """
        return (F_H2/(F_Ph0-F_Ph))

class NumericalMethods:
    
    @staticmethod
    def rk4_1step(inletPressure, inletTemperature, p, temperature, feedRate, Ac, a, Ta, F_Ph, F_H2O, F_CO, F_H2, F_CO2, phenolFraction, alpha):
        h = ReactorConstants.StepSize
        F_N2 = ReactorConstants.YI0 * feedRate
        m_calo = ReactorConstants.m_calo
        
        #k1 calculation       
        k1_temperature = temperature
        k1_p = p
        k1_Ta = Ta
        k1_F_Ph = F_Ph
        k1_F_H2O = F_H2O
        k1_F_CO = F_CO
        k1_F_H2 = F_H2
        k1_F_CO2 = F_CO2
        k1_F_T = ReactorEquations.F_T(k1_F_Ph, k1_F_H2O, k1_F_CO, k1_F_H2, k1_F_CO2, F_N2)
        
        k1_C_Ph = ReactorEquations.concentration(inletPressure, k1_p*inletPressure, inletTemperature, k1_temperature, k1_F_Ph, k1_F_T)
        k1_C_H2O = ReactorEquations.concentration(inletPressure, k1_p*inletPressure, inletTemperature, k1_temperature, k1_F_H2O, k1_F_T)
        k1_C_CO = ReactorEquations.concentration(inletPressure, k1_p*inletPressure, inletTemperature, k1_temperature, k1_F_CO, k1_F_T)
        k1_C_H2 = ReactorEquations.concentration(inletPressure, k1_p*inletPressure, inletTemperature, k1_temperature, k1_F_H2, k1_F_T)
        k1_C_CO2 = ReactorEquations.concentration(inletPressure, k1_p*inletPressure, inletTemperature, k1_temperature, k1_F_CO2, k1_F_T)
        k1_r_srp = ReactorEquations.r_srp(k1_temperature, k1_C_Ph, k1_C_H2O)
        k1_r_wgs = ReactorEquations.r_wgs(k1_temperature, k1_C_Ph, k1_C_CO, k1_C_H2O, k1_C_CO2, k1_C_H2)
        
        k1_dTdW = ReactorEquations.bilan_E(k1_temperature, k1_Ta, a, k1_F_Ph, k1_F_H2O, k1_F_H2, k1_F_CO, k1_F_CO2, F_N2, k1_C_Ph, k1_C_H2O, k1_C_CO, k1_C_CO2, k1_C_H2)
        k1_dpdW = ReactorEquations.ergun(alpha, k1_p, k1_temperature, inletTemperature, k1_F_T, feedRate)
        k1_dTadW = ReactorEquations.echange_chal(k1_temperature, k1_Ta, a, m_calo)
        k1_dF_Ph_dW = ReactorEquations.r_Ph(k1_r_srp) 
        k1_dF_H20_dW = ReactorEquations.r_H2O(k1_r_srp, k1_r_wgs) 
        k1_dF_CO_dW = ReactorEquations.r_CO(k1_r_srp, k1_r_wgs) 
        k1_dF_H2_dW = ReactorEquations.r_H2(k1_r_srp, k1_r_wgs) 
        k1_dF_CO2_dW = ReactorEquations.r_CO2(k1_r_wgs)
        
        #k2 calculation
        k2_temperature = temperature + 0.5*h*k1_dTdW
        k2_p = p + (0.5 * h * k1_dpdW)
        k2_Ta = Ta + (0.5 * h * k1_dTadW)
        k2_F_Ph = F_Ph + (h * 0.5 *k1_dF_Ph_dW)
        k2_F_H2O = F_H2O + (h * 0.5 *k1_dF_H20_dW)
        k2_F_CO = F_CO + (h * 0.5 *k1_dF_CO_dW)
        k2_F_H2 = F_H2 + (h * 0.5 *k1_dF_H2_dW)
        k2_F_CO2 = F_CO2 + (h * 0.5 *k1_dF_CO2_dW)
        k2_F_T = ReactorEquations.F_T(k2_F_Ph, k2_F_H2O, k2_F_CO, k2_F_H2, k2_F_CO2, F_N2)
        
        k2_C_Ph = ReactorEquations.concentration(inletPressure, k2_p*inletPressure, inletTemperature, k2_temperature, k2_F_Ph, k2_F_T)
        k2_C_H2O = ReactorEquations.concentration(inletPressure, k2_p*inletPressure, inletTemperature, k2_temperature, k2_F_H2O, k2_F_T)
        k2_C_CO = ReactorEquations.concentration(inletPressure, k2_p*inletPressure, inletTemperature, k2_temperature, k2_F_CO, k2_F_T)
        k2_C_H2 = ReactorEquations.concentration(inletPressure, k2_p*inletPressure, inletTemperature, k2_temperature, k2_F_H2, k2_F_T)
        k2_C_CO2 = ReactorEquations.concentration(inletPressure, k2_p*inletPressure, inletTemperature, k2_temperature, k2_F_CO2, k2_F_T)
        k2_r_srp = ReactorEquations.r_srp(k2_temperature, k2_C_Ph, k2_C_H2O)
        k2_r_wgs = ReactorEquations.r_wgs(k2_temperature, k2_C_Ph, k2_C_CO, k2_C_H2O, k2_C_CO2, k2_C_H2)

        k2_dTdW = ReactorEquations.bilan_E(k2_temperature, k2_Ta, a, k2_F_Ph, k2_F_H2O, k2_F_H2, k2_F_CO, k2_F_CO2, F_N2, k2_C_Ph, k2_C_H2O, k2_C_CO, k2_C_CO2, k2_C_H2)
        k2_dpdW = ReactorEquations.ergun(alpha, k2_p, k2_temperature, inletTemperature, k2_F_T, feedRate)
        k2_dTadW = ReactorEquations.echange_chal(k2_temperature, k2_Ta, a, m_calo)
        k2_dF_Ph_dW = ReactorEquations.r_Ph(k2_r_srp) 
        k2_dF_H2O_dW = ReactorEquations.r_H2O(k2_r_srp, k2_r_wgs) 
        k2_dF_CO_dW = ReactorEquations.r_CO(k2_r_srp, k2_r_wgs) 
        k2_dF_H2_dW = ReactorEquations.r_H2(k2_r_srp, k2_r_wgs) 
        k2_dF_CO2_dW = ReactorEquations.r_CO2(k2_r_wgs)
        
        #k3 calcualtion
        k3_temperature = temperature + (0.5 * h * k2_dTdW)
        k3_p = p + (0.5 * h * k2_dpdW)
        k3_Ta = Ta + (0.5 * h * k2_dTadW)
        k3_F_Ph = F_Ph + (0.5 * h * k2_dF_Ph_dW)
        k3_F_H2O = F_H2O + (0.5 * h * k2_dF_H2O_dW)
        k3_F_CO = F_CO + (0.5 * h * k2_dF_CO_dW)
        k3_F_H2 = F_H2 + (0.5 * h * k2_dF_H2_dW)
        k3_F_CO2 = F_CO2 + (0.5 * h * k2_dF_CO2_dW)
        k3_F_T = ReactorEquations.F_T(k3_F_Ph, k3_F_H2O, k3_F_CO, k3_F_H2, k3_F_CO2, F_N2)
        
        k3_C_Ph = ReactorEquations.concentration(inletPressure, k3_p*inletPressure, inletTemperature, k3_temperature, k3_F_Ph, k3_F_T)
        k3_C_H2O = ReactorEquations.concentration(inletPressure, k3_p*inletPressure, inletTemperature, k3_temperature, k3_F_H2O, k3_F_T)
        k3_C_CO = ReactorEquations.concentration(inletPressure, k3_p*inletPressure, inletTemperature, k3_temperature, k3_F_CO, k3_F_T)
        k3_C_H2 = ReactorEquations.concentration(inletPressure, k3_p*inletPressure, inletTemperature, k3_temperature, k3_F_H2, k3_F_T)
        k3_C_CO2 = ReactorEquations.concentration(inletPressure, k3_p*inletPressure, inletTemperature, k3_temperature, k3_F_CO2, k3_F_T)
        k3_r_srp = ReactorEquations.r_srp(k3_temperature, k3_C_Ph, k3_C_H2O)
        k3_r_wgs = ReactorEquations.r_wgs(k3_temperature, k3_C_Ph, k3_C_CO, k3_C_H2O, k3_C_CO2, k3_C_H2)

        k3_dTdW = ReactorEquations.bilan_E(k3_temperature, k3_Ta, a, k3_F_Ph, k3_F_H2O, k3_F_H2, k3_F_CO, k3_F_CO2, F_N2, k3_C_Ph, k3_C_H2O, k3_C_CO, k3_C_CO2, k3_C_H2)
        k3_dpdW = ReactorEquations.ergun(alpha, k3_p, k3_temperature, inletTemperature, k3_F_T, feedRate)
        k3_dTadW = ReactorEquations.echange_chal(k3_temperature, k3_Ta, a, m_calo)
        k3_dF_Ph_dW = ReactorEquations.r_Ph(k3_r_srp) 
        k3_dF_H2O_dW = ReactorEquations.r_H2O(k3_r_srp, k3_r_wgs) 
        k3_dF_CO_dW = ReactorEquations.r_CO(k3_r_srp, k3_r_wgs) 
        k3_dF_H2_dW = ReactorEquations.r_H2(k3_r_srp, k3_r_wgs) 
        k3_dF_CO2_dW = ReactorEquations.r_CO2(k3_r_wgs)

        #k4 calculation
        k4_temperature = temperature + (h * k3_temperature)
        k4_p = p + (h * k3_p)
        k4_Ta = Ta + (h * k3_Ta)
        k4_F_Ph = F_Ph + (h * k3_F_Ph)
        k4_F_H2O = F_H2O + (h * k3_F_H2O)
        k4_F_CO = F_CO + (h * k3_F_CO)
        k4_F_H2 = F_H2 + (h * k3_F_H2)
        k4_F_CO2 = F_CO2 + (h * k3_F_CO2)
        k4_F_T = ReactorEquations.F_T(k4_F_Ph, k4_F_H2O, k4_F_CO, k4_F_H2, k4_F_CO2, F_N2)
        
        k4_C_Ph = ReactorEquations.concentration(inletPressure, k4_p*inletPressure, inletTemperature, k4_temperature, k4_F_Ph, k4_F_T)
        k4_C_H2O = ReactorEquations.concentration(inletPressure, k4_p*inletPressure, inletTemperature, k4_temperature, k4_F_H2O, k4_F_T)
        k4_C_CO = ReactorEquations.concentration(inletPressure, k4_p*inletPressure, inletTemperature, k4_temperature, k4_F_CO, k4_F_T)
        k4_C_H2 = ReactorEquations.concentration(inletPressure, k4_p*inletPressure, inletTemperature, k4_temperature, k4_F_H2, k4_F_T)
        k4_C_CO2 = ReactorEquations.concentration(inletPressure, k4_p*inletPressure, inletTemperature, k4_temperature, k4_F_CO2, k4_F_T)
        k4_r_srp = ReactorEquations.r_srp(k4_temperature, k4_C_Ph, k4_C_H2O)
        k4_r_wgs = ReactorEquations.r_wgs(k4_temperature, k4_C_Ph, k4_C_CO, k4_C_H2O, k4_C_CO2, k4_C_H2)

        k4_dTdW = ReactorEquations.bilan_E(k4_temperature, k4_Ta, a, k4_F_Ph, k4_F_H2O, k4_F_H2, k4_F_CO, k4_F_CO2, F_N2, k4_C_Ph, k4_C_H2O, k4_C_CO, k4_C_CO2, k4_C_H2)
        k4_dpdW = ReactorEquations.ergun(alpha, k4_p, k4_temperature, inletTemperature, k4_F_T, feedRate)
        k4_dTadW = ReactorEquations.echange_chal(k4_temperature, k4_Ta, a, m_calo)
        k4_dF_Ph_dW = ReactorEquations.r_Ph(k4_r_srp) 
        k4_dF_H2O_dW = ReactorEquations.r_H2O(k4_r_srp, k4_r_wgs) 
        k4_dF_CO_dW = ReactorEquations.r_CO(k4_r_srp, k4_r_wgs) 
        k4_dF_H2_dW = ReactorEquations.r_H2(k4_r_srp, k4_r_wgs) 
        k4_dF_CO2_dW = ReactorEquations.r_CO2(k4_r_wgs)
        
        #calculation of final ROC
        dTdW = (1.0/6.0)*(k1_dTdW + 2*k2_dTdW + 2*k3_dTdW + k4_dTdW)
        dpdW = (1.0/6.0)*(k1_dpdW + 2*k2_dpdW + 2*k3_dpdW + k4_dpdW)
        dTadW = (1.0/6.0)*(k1_dTadW + 2*k2_dTadW + 2*k3_dTadW + k4_dTadW)
        dF_Ph_dW = (1.0/6.0)*(k1_dF_Ph_dW + 2*k2_dF_Ph_dW + 2*k3_dF_Ph_dW + k4_dF_Ph_dW)
        dF_H2O_dW = (1.0/6.0)*(k1_dF_H20_dW + 2*k2_dF_H2O_dW + 2*k3_dF_H2O_dW + k4_dF_H2O_dW)
        dF_CO_dW = (1.0/6.0)*(k1_dF_CO_dW + 2*k2_dF_CO_dW + 2*k3_dF_CO_dW + k4_dF_CO_dW)
        dF_H2_dW = (1.0/6.0)*(k1_dF_H2_dW + 2*k2_dF_H2_dW + 2*k3_dF_H2_dW + k4_dF_H2_dW)
        dF_CO2_dW = (1.0/6.0)*(k1_dF_CO2_dW + 2*k2_dF_CO2_dW + 2*k3_dF_CO2_dW + k4_dF_CO2_dW)
        
        return dTdW, dpdW, dTadW, dF_Ph_dW, dF_H2O_dW, dF_CO_dW, dF_H2_dW, dF_CO2_dW

def GraphRun(t, pres, fr, pf, ar, ss):
    inletTemperature =  t#K
    inletPressure = pres #atm
    feedRate = fr #mol/s
    phenolFraction = pf #molPh/mol
    Ac = ar #m^2

    h = ss #Step Size

    ####################################################################################################################################################################
    #Constant and Array Declaration
    ####################################################################################################################################################################
    a = ReactorEquations.a(Ac)
    temperature = inletTemperature
    p = 1
    Ta = ReactorConstants.T_calo_in
    F_Ph = feedRate * phenolFraction
    F_H2O = feedRate * (1-phenolFraction-ReactorConstants.YI0)
    F_CO = 0
    F_H2 = 0
    F_CO2 = 0
    F_N2 = feedRate * ReactorConstants.YI0
    F_T = ReactorEquations.F_T(F_Ph, F_H2O, F_CO, F_H2, F_CO2, F_N2)
    alpha = ReactorEquations.alpha(inletPressure, inletTemperature, temperature, phenolFraction, Ac, inletPressure, F_T, phenolFraction)

    C_Ph_arr = [ReactorEquations.concentration(inletPressure,inletPressure, inletTemperature, inletTemperature, F_Ph, F_T)]
    C_H2O_arr = [ReactorEquations.concentration(inletPressure,inletPressure, inletTemperature, inletTemperature, F_H2O, F_T)]
    C_CO_arr = [ReactorEquations.concentration(inletPressure,inletPressure, inletTemperature, inletTemperature, F_CO, F_T)]
    C_H2_arr = [ReactorEquations.concentration(inletPressure,inletPressure, inletTemperature, inletTemperature, F_H2, F_T)]
    C_CO2_arr = [ReactorEquations.concentration(inletPressure,inletPressure, inletTemperature, inletTemperature, F_CO2, F_T)]
    C_N2_arr = [ReactorEquations.concentration(inletPressure,inletPressure, inletTemperature, inletTemperature, F_N2, F_T)]

    F_Ph_arr = [F_Ph]
    F_H2O_arr = [F_H2O]
    F_CO_arr = [F_CO]
    F_H2_arr = [F_H2]
    F_CO2_arr = [F_CO2]
    F_N2_arr = [F_N2]

    catalystWeigtharr = [0]
    conversionarr = [0]

    ####################################################################################################################################################################
    #RK4 implementation usage
    ####################################################################################################################################################################
    catalystWeigth = 0
    while(F_H2 < 25.0 and catalystWeigth < 30 and p>0):
        
        TdW, dpdW, dTadW, dF_Ph_dW, dF_H2O_dW, dF_CO_dW, dF_H2_dW, dF_CO2_dW = \
        NumericalMethods.rk4_1step(inletPressure, inletTemperature, p, temperature, feedRate, Ac, a, Ta, F_Ph, F_H2O, F_CO, F_H2, F_CO2, phenolFraction, alpha)
        
        temperature += h * TdW
        p += h * dpdW
        Ta += h * dTadW
        F_Ph += h * dF_Ph_dW
        F_H2O += h * dF_H2O_dW
        F_CO += h * dF_CO_dW
        F_H2 += h * dF_H2_dW
        F_CO2 += h * dF_CO2_dW
        F_T = ReactorEquations.F_T(F_Ph, F_H2O, F_CO, F_H2, F_CO2, F_N2)
        
        C_Ph_arr.append(ReactorEquations.concentration(inletPressure, p*inletPressure, inletTemperature, temperature, F_Ph, F_T))
        C_H2O_arr.append(ReactorEquations.concentration(inletPressure, p*inletPressure, inletTemperature, temperature, F_H2O, F_T))
        C_CO_arr.append(ReactorEquations.concentration(inletPressure, p*inletPressure, inletTemperature, temperature, F_CO, F_T))
        C_H2_arr.append(ReactorEquations.concentration(inletPressure, p*inletPressure, inletTemperature, temperature, F_H2, F_T))
        C_CO2_arr.append(ReactorEquations.concentration(inletPressure, p*inletPressure, inletTemperature, temperature, F_CO2, F_T))
        C_N2_arr.append(ReactorEquations.concentration(inletPressure, p*inletPressure, inletTemperature, temperature, F_N2, F_T))
        
        F_Ph_arr.append(F_Ph)
        F_H2O_arr.append(F_H2O)
        F_CO_arr.append(F_CO)
        F_H2_arr.append(F_H2)
        F_CO2_arr.append(F_CO2)
        F_N2_arr.append(F_N2)

        catalystWeigth += h
        catalystWeigtharr.append(catalystWeigth)
        conversionarr.append(ReactorEquations.conversion(feedRate*phenolFraction, F_Ph))


    ####################################################################################################################################################################
    #Error Handling and Plotting
    ####################################################################################################################################################################
    if(F_H2 < 25.0):
        print("Non Convering Solution")
    else:
        concFig, concAxis = plt.subplots(1, 1)
        concFig.set_figheight(10)
        concFig.set_figwidth(10)
        concFig.tight_layout(pad=8.0)

        concAxis.plot(catalystWeigtharr, C_Ph_arr, c='blue', label='Phenol')
        concAxis.plot(catalystWeigtharr, C_H2O_arr, c='orange', label='H2O')
        concAxis.plot(catalystWeigtharr, C_CO_arr, c='red', label='CO')
        concAxis.plot(catalystWeigtharr, C_H2_arr, c='purple', label='H2')
        concAxis.plot(catalystWeigtharr, C_CO2_arr, c='cyan', label='CO2')
        concAxis.plot(catalystWeigtharr, C_N2_arr, c='pink', label='N2')
        concAxis.legend()
        concAxis.set_title("Concentration des especes en fonction de la masse du catalyseur")
        concAxis.set_xlabel('Masse du Catalyseur [kg]')
        concAxis.set_ylabel('Concentration [mol/m^3]')


        flowFig, flowAxis = plt.subplots(1, 1)
        flowFig.set_figheight(10)
        flowFig.set_figwidth(10)
        flowFig.tight_layout(pad=8.0)

        flowAxis.plot(catalystWeigtharr, F_Ph_arr, c='blue', label='Phenol')
        flowAxis.plot(catalystWeigtharr, F_H2O_arr, c='orange', label='H2O')
        flowAxis.plot(catalystWeigtharr, F_CO_arr, c='red', label='CO')
        flowAxis.plot(catalystWeigtharr, F_H2_arr, c='purple', label='H2')
        flowAxis.plot(catalystWeigtharr, F_CO2_arr, c='cyan', label='CO2')
        flowAxis.plot(catalystWeigtharr, F_N2_arr, c='pink', label='N2')
        flowAxis.legend()
        flowAxis.set_title("Debit molaire des especes en fonction de la masse du catalyseur")
        flowAxis.set_xlabel('Masse du Catalyseur [kg]')
        flowAxis.set_ylabel('Debit Molaire [mol/s]')

        convFig, convAxis = plt.subplots(1, 1)
        convFig.set_figheight(10)
        convFig.set_figwidth(10)
        convFig.tight_layout(pad=8.0)

        convAxis.plot(catalystWeigtharr, conversionarr, c='purple', label='Conversion')
        convAxis.legend()
        convAxis.set_title("Conversion de phenol en fonction de la masse du catalyseur")
        convAxis.set_xlabel('Masse du Catalyseur [kg]')
        convAxis.set_ylabel('Conversion')

        plt.show()