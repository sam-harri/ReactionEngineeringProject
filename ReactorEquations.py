import matplotlib.pyplot as plt
import numpy as np
from ReactorConstants import ReactorConstants
from typing import List
from math import sqrt

FA0: float = ReactorConstants.FA0
alpha: float = ReactorConstants.alpha
weigth: float = ReactorConstants.weight

class ReactorEquations:
    
    @staticmethod
    def pressureDropAnalytical():
        return sqrt(1-alpha*weigth)
    
    @staticmethod
    def levenspielPlot():
        xdata : List[float] = np.linspace(0,1)
        ydata : List[float] = [FA0/ReactorEquations.calcDissapearence(conversion) for conversion in xdata]
        plt.plot(ydata,xdata)
        plt.xlabel("Conversion")
        plt.ylabel("FA0/-Ra")
        plt.show()
        
    @staticmethod
    def calcDissapearence(conversion : float):
        pass