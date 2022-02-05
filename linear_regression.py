import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv("data.csv")
#data.plot()
#data.plot(kind='density')
X = data.iloc[:, 0]
Y = data.iloc[:, 1]
#plt.scatter(X,Y)
#plt.show()
m = 0
c = 0
L = 0.0001
epochs = 1000
n = float(len(X))
for i in range(epochs):
    Y_pred = m*X + c
    D_m = (-2/n) * sum(X* (Y - Y_pred))
    D_c = (-2/n)*sum(Y - Y_pred)
    m = m - L*D_m
    c = c - L*D_c
    
print (m, c)


plt.scatter(X,Y)
Y_pred = m*X + c
plt.plot([min(X),max(X)],[min(Y_pred), max(Y_pred)], color='red')
plt.show()
