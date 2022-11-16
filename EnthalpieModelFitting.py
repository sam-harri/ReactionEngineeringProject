from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

x_surf, y_surf = np.meshgrid(
    np.linspace(1,10,100),
    np.linspace(600,800,100)
)
surfaceX = pd.DataFrame({'pressure' : x_surf.ravel(), 'temperature' : y_surf.ravel()})

df = pd.DataFrame([
    [1,600,8571.0],
    [1,610,8716.4],
    [1,620,8861.9],
    [1,630,9007.4],
    [1,640,9153.0],
    [1,650,9298.6],
    [1,660,9444.3],
    [1,670,9590.0],
    [1,680,9735.8],
    [1,690,9881.7],
    [1,700,10027.7],
    [1,710,10173.8],
    [1,720,10320.0],
    [1,730,10466.2],
    [1,740,10612.6],
    [1,750,10759.1],
    [1,760,10905.7],
    [1,770,11052.3],
    [1,780,11199.1],
    [1,790,11346.1],
    [1,800,11493.1],
    [5,600,8573.6],
    [5,610,8719.1],
    [5,620,8864.6],
    [5,630,9010.1],
    [5,640,9155.7],
    [5,650,9301.3],
    [5,660,9447.0],
    [5,670,9592.8],
    [5,680,9738.6],
    [5,690,9884.6],
    [5,700,10030.6],
    [5,710,10176.7],
    [5,720,10322.9],
    [5,730,10469.1],
    [5,740,10615.5],
    [5,750,10762.0],
    [5,760,10908.6],
    [5,770,11055.3],
    [5,780,11202.1],
    [5,790,11349.1],
    [5,800,11496.1],
    [10,600,8576.6],
    [10,610,8722.4],
    [10,620,8868.0],
    [10,630,9013.5],
    [10,640,9159.1],
    [10,650,9304.8],
    [10,660,9450.5],
    [10,670,9596.3],
    [10,680,9742.2],
    [10,690,9888.1],
    [10,700,10034.2],
    [10,710,10180.3],
    [10,720,10326.5],
    [10,730,10472.8],
    [10,740,10619.2],
    [10,750,10765.7],
    [10,760,10912.3],
    [10,770,11059.0],
    [10,780,11205.9],
    [10,790,11352.8],
    [10,800,11499.9]
    ], columns=['pressure', 'temperature', 'enthalpy'])

model = make_pipeline(
    PolynomialFeatures(degree=2),
    LinearRegression()
)

model.fit(df[['pressure', 'temperature']], df[['enthalpy']])

predictedSurface = model.predict(surfaceX)

figure = plt.figure(figsize=(15,10))

axis = figure.add_subplot(221, projection='3d')
axis.scatter(df['pressure'],df['temperature'], df['enthalpy'], c='mediumorchid', marker='o', alpha=0.5)
axis.plot_surface(x_surf, y_surf, predictedSurface.reshape(x_surf.shape), alpha=0.3, color='mediumorchid')
axis.set_title("Multiple Polynomal Regression Model Hydrogen Enthalpy")
axis.set_xlabel("Pressure [atm]")
axis.set_ylabel("Temperature [K]")
axis.set_zlabel("Enthalpy [J/g]")

axis2 = figure.add_subplot(222, projection='3d')


plt.show()