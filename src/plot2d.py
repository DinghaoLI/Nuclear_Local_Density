import scipy.io as sio  
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

X = np.loadtxt("rVals.txt")
Y = np.loadtxt("zVals.txt")
Z = np.loadtxt("plot2d.txt")

fig = plt.figure()
surf1 = plt.contourf(X, Y, Z, 1000)
fig.colorbar(surf1)
plt.xlabel('R')
plt.ylabel('Z')
plt.title('Density')
plt.show()