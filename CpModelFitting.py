from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np

base = np.linspace(800,1100)

#-0.00023050000x^2;0.51767215x^1;-268.3843301362517x^0
xPhenol = np.array([800, 900, 1000,1100])
yPhenol = np.array([212.14, 223.19, 232.49, 240.41])

#kSRP model training
poly_PhenolCp = PolynomialFeatures(degree=2, include_bias=False)
poly_features_PhenolCp = poly_PhenolCp.fit_transform(xPhenol.reshape(-1, 1))
poly_reg_model_PhenolCp = LinearRegression()
poly_reg_model_PhenolCp.fit(poly_features_PhenolCp, yPhenol)

#Running data through predicted model to plot polynomial
ySRP_train = poly_reg_model_PhenolCp.predict(poly_features_PhenolCp)
ySRP_test = poly_reg_model_PhenolCp.predict(poly_PhenolCp.fit_transform(base.reshape(-1, 1)))

figure, axis = plt.subplots(1, 1)
figure.set_figheight(10)
figure.set_figwidth(10)
figure.tight_layout(pad=8.0)

axis.scatter(xPhenol,yPhenol, color='mediumorchid')

axis.plot(base, ySRP_test, c='mediumorchid', label='SRP')
axis.legend()
axis.set_title("LSS Polynomial Regression Models of Cp of Phenol vs T")
axis.set_xlabel('Temperature (K)')
axis.set_ylabel('Cp (J/mol*K)')
plt.show()

# print(poly_reg_model_PhenolCp.coef_)
# print(poly_reg_model_PhenolCp.intercept_)