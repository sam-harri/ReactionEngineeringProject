import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from ReactorConstants import ReactorConstants
from ReactorEquations import ReactorEquations
from NumericalMethods import NumericalMethods
import matplotlib.pyplot as plt

####################################################################################################################################################################
#Run Parameters
####################################################################################################################################################################
inletTemperature = 1073.15 #K
inletPressure = 15 #atm
feedRate = 150 #mol/s
phenolFraction = 0.0325 #molPh/mol
Ac = 0.02 #m^2

h = 0.0001 #Step Size

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