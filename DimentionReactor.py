from ReactorConstants import ReactorConstants
from ReactorEquations import ReactorEquations
from NumericalMethods import NumericalMethods

class DimentionReactor:
    
    @staticmethod
    def dimentionlizeReactor(inletTemperature: float, inletPressure: float, feedRate: float, phenolFraction: float, Ac: float):
        
        h = ReactorConstants.StepSize
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
        
        catalystWeigth = 0
        while(F_H2 < 25.0 and catalystWeigth < 1000):
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
            
            catalystWeigth += h
        
        conversion = ReactorEquations.conversion(feedRate*phenolFraction, F_Ph)
        gloSelect = ReactorEquations.select_glo(F_H2, F_CO)
        
        if(F_H2 > 25.0):
            return catalystWeigth, conversion, gloSelect
        return 0,0,0