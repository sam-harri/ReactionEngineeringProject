from typing import Callable
from ReactorEquations import ReactorEquations
from ReactorConstants import ReactorConstants

class NumericalMethods:
    
    @staticmethod
    def bisection(upper: float, lower: float, func: Callable):
        """
        Returns the root of func found between [Upper, Lower]
        """
        old: float = 1 #temp set variables, get overridden after first loop
        new: float = 0
        middle: float = (upper + lower)/2
        while(abs(old-new)>0.00001):
            if(func(middle)*func(upper)>0):
                upper = middle
            else:
                lower = middle
            old = func(middle)
            middle = (upper + lower)/2
            new = func(middle)
        return middle

    @staticmethod
    def rungeKutta(x0: float, y0: float, x: float, h: float, dydx: Callable):
        """
        Evaluates ODE dydx at point x with a step size of h
        using the initial conditions x0, y0
        """
        n: int = (int)((x - x0)/h)
        y: float = y0
        for _ in range(1, n + 1):
            k1: float = h * dydx(y)
            k2: float = h * dydx(y + 0.5 * k1)
            k3: float = h * dydx(y + 0.5 * k2)
            k4: float = h * dydx(y + k3)
            y +=(1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
            x0 += h
        return y
    
    def rk4_dpdW(alpha, p0, temperature, inletTemperature, F_T, feedrate):
        h = ReactorConstants.StepSize
        k1 = h * ReactorEquations.ergun(alpha, p0, temperature, inletTemperature, F_T, feedrate)
        k2 = h * ReactorEquations.ergun(alpha, p0 + 0.5*k1, temperature, inletTemperature, F_T, feedrate)
        k3 = h * ReactorEquations.ergun(alpha, p0 + 0.5*k2, temperature, inletTemperature, F_T, feedrate)
        k4 = h * ReactorEquations.ergun(alpha, p0 + k3, temperature, inletTemperature, F_T, feedrate)
        return (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
    
    def rk4_dTdW(temperature, Ta, a, F_Ph, F_H2O, F_H2, F_CO, F_CO2, F_N2, C_Ph, C_H2O, C_CO, C_CO2, C_H2):
        h = ReactorConstants.StepSize
        k1 = h * ReactorEquations.bilan_E(temperature, Ta, a, F_Ph, F_H2O, F_H2, F_CO, F_CO2, F_N2, C_Ph, C_H2O, C_CO, C_CO2, C_H2)
        k2 = h * ReactorEquations.bilan_E(temperature + 0.5*k1, Ta, a, F_Ph, F_H2O, F_H2, F_CO, F_CO2, F_N2, C_Ph, C_H2O, C_CO, C_CO2, C_H2)
        k3 = h * ReactorEquations.bilan_E(temperature + 0.5*k2, Ta, a, F_Ph, F_H2O, F_H2, F_CO, F_CO2, F_N2, C_Ph, C_H2O, C_CO, C_CO2, C_H2)
        k4 = h * ReactorEquations.bilan_E(temperature + k3, Ta, a, F_Ph, F_H2O, F_H2, F_CO, F_CO2, F_N2, C_Ph, C_H2O, C_CO, C_CO2, C_H2)
        return (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
    
    def rk4_dTadW(temperature, Ta, a, m_calo):
        h = ReactorConstants.StepSize
        k1 = h * ReactorEquations.echange_chal(temperature, Ta, a, m_calo)
        k2 = h * ReactorEquations.echange_chal(temperature, Ta + 0.5*k1, a, m_calo)
        k3 = h * ReactorEquations.echange_chal(temperature, Ta + 0.5*k2, a, m_calo)
        k4 = h * ReactorEquations.echange_chal(temperature, Ta + k3, a, m_calo)
        return (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)