from Helper.Polynomial import Polynomial

class ReactorConstants:
    
    #add all constants here so that we can easily update them acrossall classes
    #changes here will be reflected on everything else
    #make sure to delcare these as the constants by using the class name everywhere else though
    #ex : FA0: float = ReactorConstants.FA0
    
    FA0: float = 25 #molH2/s à la sortie
    YI0: float = 0.1 #molI/mol
    U: float = 16.0 #W/m**2K
    alpha: float = 0 #temp
    weight : float = 0 #temp
    KP: Polynomial = Polynomial("0.000041x^2;-0.0600983x^1;22.88680282250014x^0")
    KH2O: Polynomial = Polynomial("0.0000745x^2;-0.11574935x^1;50.048421226250284x^0")
    Deq: float = 0.005 #m diamètre équivalent catalyseur
    epsilon: float = 0.39 #porosité
    rho: float = 5740 #kg/m^3 masse volumique catalyseur