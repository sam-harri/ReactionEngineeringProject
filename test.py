from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np

# xSRP = np.array([873.15, 973.15, 1073.15])
# ySRP = np.array([7.89, 17.10, 21.70])

x = np.array([873.15, 973.15, 1073.15])
y = np.array([7.89, 17.10, 21.70])

xKGS = np.array([873.15, 973.15, 1073.15])
yKGS = np.array([6.45, 13.90, 17.81])

xRWGS = [873.15, 973.15, 1073.15]
yRWGS = [5.10, 10.63, 13.65]

poly = PolynomialFeatures(degree=2, include_bias=False)
poly_features = poly.fit_transform(x.reshape(-1, 1))
poly_reg_model = LinearRegression()
poly_reg_model.fit(poly_features, y)

y_predicted = poly_reg_model.predict(poly_features)
y_test = poly_reg_model.predict(poly.fit_transform(np.linspace(800,1100).reshape(-1, 1)))

# depicting the polynomial graph
plt.figure(figsize=(10, 6))
plt.scatter(x,y)
plt.plot(np.linspace(800,1100), y_test, c='red')
plt.show()

