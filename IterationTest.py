from ReactorConstants import ReactorConstants
from ReactorEquations import ReactorEquations
from NumericalMethods import NumericalMethods
import csv
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

h = 0.001
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

with open('testCsv.csv','w', newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["catalystWeigth", "temperature", "p", "Ta", "F_Ph", "F_H2O", "F_CO", "F_H2", "F_CO2"])

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
    
    csvArr = [catalystWeigth, temperature, p, Ta, F_Ph, F_H2O, F_CO, F_H2, F_CO2]
    with open('testCsv.csv','a', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(csvArr)
    
    catalystWeigth += h
    
end = time.time()

if(F_H2 < 25.0):
    print("poop")

print("The time of execution of above program is :",
      (end-start) * 10**3, "ms")