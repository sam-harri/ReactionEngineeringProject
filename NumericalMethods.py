from typing import Callable

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
            k1: float = h * dydx(x0, y)
            k2: float = h * dydx(x0 + 0.5 * h, y + 0.5 * k1)
            k3: float = h * dydx(x0 + 0.5 * h, y + 0.5 * k2)
            k4: float = h * dydx(x0 + h, y + k3)
            y +=(1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
            x0 += h
        return y