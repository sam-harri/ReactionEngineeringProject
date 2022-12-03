from ReactorEquations import ReactorEquations
from ReactorConstants import ReactorConstants

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
        k4_temperature = temperature + (h * 0.5 * k3_temperature)
        k4_p = p + (h * 0.5 * k3_p)
        k4_Ta = Ta + (h * 0.5 * k3_Ta)
        k4_F_Ph = F_Ph + (h * 0.5 * k3_F_Ph)
        k4_F_H2O = F_H2O + (h * 0.5 * k3_F_H2O)
        k4_F_CO = F_CO + (h * 0.5 * k3_F_CO)
        k4_F_H2 = F_H2 + (h * 0.5 * k3_F_H2)
        k4_F_CO2 = F_CO2 + (h * 0.5 * k3_F_CO2)
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