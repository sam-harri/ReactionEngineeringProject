
class Term:
    def __init__(self,coefficient: float, power: int) -> None:
        self.coefficient = coefficient
        self.power = power
        
    def evaluate(self,xval: float):
        return self.coefficient * (pow(xval,self.power))