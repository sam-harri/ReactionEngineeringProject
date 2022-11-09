from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np

# xSRP = [873.15, 973.15, 1073.15]
# ySRP = [7.89, 17.10, 21.70]

x = [873.15, 973.15, 1073.15]
y = [7.89, 17.10, 21.70]

xKGS = [873.15, 973.15, 1073.15]
yKGS = [6.45, 13.90, 17.81]

xRWGS = [873.15, 973.15, 1073.15]
yRWGS = [5.10, 10.63, 13.65]

poly = PolynomialFeatures(degree=2, include_bias=False)
# creating a new feature
poly_features = poly.fit_transform(x.reshape(-1, 1))
# creating a polynomial regression model
from sklearn.linear_model import LinearRegression
poly_reg_model = LinearRegression()
poly_reg_model.fit(poly_features, y)
y_predicted = poly_reg_model.predict(poly_features)
# depicting the polynomial graph
plt.figure(figsize=(10, 6))
plt.scatter(x,y)
plt.plot(x, y_predicted, c='red')
plt.show()

