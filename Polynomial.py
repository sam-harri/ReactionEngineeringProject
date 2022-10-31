from typing import List
from Term import Term


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
    
    def integrate(self):
        for term in self.termArr:
            term.coefficient /= term.power+1
            term.power += 1