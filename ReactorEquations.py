import matplotlib.pyplot as plt
import numpy as np
from ReactorConstants import ReactorConstants
from typing import List

class ReactorEquations:
    
    @staticmethod
    def levenspielPlot():
        ydata : List[float] = np.linspace(0,1)
        xdata : List[float] = [ReactorEquations.calcDissapearence(conversion) for conversion in ydata]
        plt.plot(ydata,xdata)
        plt.show()
    
    @staticmethod
    def calcDissapearence(conversion : float):
        pass