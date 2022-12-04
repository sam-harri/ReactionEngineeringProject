import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("testCsv.csv")


figure, axis = plt.subplots(1, 1)
figure.set_figheight(10)
figure.set_figwidth(10)
figure.tight_layout(pad=8.0)

axis.scatter(df["catalystWeigth"], df["F_H2"], color='mediumorchid')

axis.legend()
axis.set_title("poop")
axis.set_xlabel('Catalyst Weigth [kg]')
axis.set_ylabel('F_H2 [mol/s]')
plt.show()