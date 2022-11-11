from Helper.Polynomial import Polynomial

class ReactorConstants:
    
    #add all constants here so that we can easily update them acrossall classes
    #changes here will be reflected on everything else
    #make sure to delcare these as the constants by using the class name everywhere else though
    #ex : FA0: float = ReactorConstants.FA0
    
    FA0: float = 25 #molA/s
    YI0: float = 0.1 #molI/mol
    U: float = 16.0 #W/m**2K
    alpha: float = 0 #temp
    weight : float = 0 #temp
    kSRP: Polynomial = Polynomial("-0.0002305000x^2;0.51767215x^1;-268.3843301362517x^0")
    kWGS: Polynomial = Polynomial("-0.0001770000x^2;0.40129510x^1;-208.9976232825012x^0")
    kRWGS: Polynomial = Polynomial("-0.00012550000x^2;0.28701065x^1;149.82328827375082x^0")