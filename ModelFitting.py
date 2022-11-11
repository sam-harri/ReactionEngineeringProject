from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np

base = np.linspace(800,1100)

#-0.00023050000x^2;0.51767215x^1;-268.3843301362517x^0
xSRP = np.array([873.15, 973.15, 1073.15])
ySRP = np.array([7.89, 17.10, 21.70])

#-0.0001770000x^2;0.4012951x^1;-208.9976232825012x^0
xWGS = np.array([873.15, 973.15, 1073.15])
yWGS = np.array([6.45, 13.90, 17.81])

#-0.00012550000x^2;0.28701065x^1;149.82328827375082x^0
xRWGS = np.array([873.15, 973.15, 1073.15])
yRWGS = np.array([5.10, 10.63, 13.65])

#kSRP model training
poly_SRP = PolynomialFeatures(degree=2, include_bias=False)
poly_features_SRP = poly_SRP.fit_transform(xSRP.reshape(-1, 1))
poly_reg_model_SRP = LinearRegression()
poly_reg_model_SRP.fit(poly_features_SRP, ySRP)

#kWGS model training
poly_WGS = PolynomialFeatures(degree=2, include_bias=False)
poly_features_WGS = poly_WGS.fit_transform(xWGS.reshape(-1, 1))
poly_reg_model_WGS = LinearRegression()
poly_reg_model_WGS.fit(poly_features_WGS, yWGS)

#kRWGS model training
poly_RWGS = PolynomialFeatures(degree=2, include_bias=False)
poly_features_RWGS = poly_RWGS.fit_transform(xRWGS.reshape(-1, 1))
poly_reg_model_RWGS = LinearRegression()
poly_reg_model_RWGS.fit(poly_features_RWGS, yRWGS)


#Running data through predicted model to plot polynomial
ySRP_train = poly_reg_model_SRP.predict(poly_features_SRP)
ySRP_test = poly_reg_model_SRP.predict(poly_SRP.fit_transform(base.reshape(-1, 1)))

yWGS_train = poly_reg_model_WGS.predict(poly_features_WGS)
yWGS_test = poly_reg_model_WGS.predict(poly_WGS.fit_transform(base.reshape(-1, 1)))

yRWGS_train = poly_reg_model_RWGS.predict(poly_features_RWGS)
yRWGS_test = poly_reg_model_RWGS.predict(poly_RWGS.fit_transform(base.reshape(-1, 1)))

plt.figure(figsize=(11, 7))

plt.scatter(xSRP,ySRP, color='mediumorchid')
plt.scatter(xWGS, yWGS, color='mediumslateblue')
plt.scatter(xRWGS, yRWGS, color='lightseagreen')

plt.plot(base, ySRP_test, c='mediumorchid', label='SRP')
plt.plot(base, yWGS_test, c='mediumslateblue', label='WGS')
plt.plot(base, yRWGS_test, c='lightseagreen', label='RWGS')
plt.legend()

plt.title("LSS Polynomial Regression Models of Kinetic Parameters of SRP Reactions using NTZ Catalyst")
plt.xlabel('Temperature (K)')
plt.ylabel('Rate Constant (mol / g-cat*h)')

plt.show()

#test point, copy paste code below and change point for model prediction
point = 873.15
SRP_estimation = poly_reg_model_SRP.predict(poly_SRP.fit_transform([[point]]))
WGS_estimation = poly_reg_model_WGS.predict(poly_WGS.fit_transform([[point]]))
RWGS_estimation = poly_reg_model_RWGS.predict(poly_RWGS.fit_transform([[point]]))


# print(poly_reg_model_RWGS.coef_)
# print(poly_reg_model_RWGS.intercept_)

