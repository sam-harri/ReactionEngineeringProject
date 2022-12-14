from ReactorConstants import ReactorConstants
from ReactorEquations import ReactorEquations
from NumericalMethods import NumericalMethods
import matplotlib.pyplot as plt
import time

start = time.time()

inletTemperature = 1023.15 #K
inletPressure = 10 #atm
feedRate = 150 #mol/s
phenolFraction = 0.04 #molPh/mol
Ac = 0.05 #m^2

# inletTemperature = 973.15 #K
# inletPressure = 3 #atm
# feedRate = 175 #mol/s
# phenolFraction = 0.02 #molPh/mol
# Ac = 0.4 #m^2

h = 0.0001
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

# with open('testCsv.csv','w', newline="") as f:
#     writer = csv.writer(f)
#     writer.writerow(["catalystWeigth", "temperature", "p", "Ta", "F_Ph", "F_H2O", "F_CO", "F_H2", "F_CO2"])

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
    
    # csvArr = [catalystWeigth, temperature, p, Ta, F_Ph, F_H2O, F_CO, F_H2, F_CO2]
    # with open('testCsv.csv','a', newline="") as f:
    #     writer = csv.writer(f)
    #     writer.writerow(csvArr)
    
    catalystWeigth += h
    catalystWeigtharr.append(catalystWeigth)
    conversionarr.append(ReactorEquations.conversion(feedRate*phenolFraction, F_Ph))
    
end = time.time()

if(F_H2 < 25.0):
    print("poop")

print("The time of execution of above program is :",
      (end-start) * 10**3, "ms")


figure, axis = plt.subplots(3, 1)
figure.set_figheight(10)
figure.set_figwidth(10)
figure.tight_layout(pad=8.0)

axis[0].plot(catalystWeigtharr, C_Ph_arr, c='blue', label='Phenol')
axis[0].plot(catalystWeigtharr, C_H2O_arr, c='orange', label='H2O')
axis[0].plot(catalystWeigtharr, C_CO_arr, c='red', label='CO')
axis[0].plot(catalystWeigtharr, C_H2_arr, c='purple', label='H2')
axis[0].plot(catalystWeigtharr, C_CO2_arr, c='cyan', label='CO2')
axis[0].plot(catalystWeigtharr, C_N2_arr, c='pink', label='N2')
axis[0].legend()
axis[0].set_title("Concentration des especes en fonction de la masse du catalyseur")
axis[0].set_xlabel('Masse du Catalyseur [kg]')
axis[0].set_ylabel('Concentration [mol/m^3]')


axis[1].plot(catalystWeigtharr, F_Ph_arr, c='blue', label='Phenol')
axis[1].plot(catalystWeigtharr, F_H2O_arr, c='orange', label='H2O')
axis[1].plot(catalystWeigtharr, F_CO_arr, c='red', label='CO')
axis[1].plot(catalystWeigtharr, F_H2_arr, c='purple', label='H2')
axis[1].plot(catalystWeigtharr, F_CO2_arr, c='cyan', label='CO2')
axis[1].plot(catalystWeigtharr, F_N2_arr, c='pink', label='N2')
axis[1].legend()
axis[1].set_title("Debit molaire des especes en fonction de la masse du catalyseur")
axis[1].set_xlabel('Masse du Catalyseur [kg]')
axis[1].set_ylabel('Debit Molaire [mol/s]')

axis[2].plot(catalystWeigtharr, conversionarr, c='purple', label='Conversion')
axis[2].legend()
axis[2].set_title("Conversion de phenol en fonction de la masse du catalyseur")
axis[2].set_xlabel('Masse du Catalyseur [kg]')
axis[2].set_ylabel('Conversion')

plt.show()