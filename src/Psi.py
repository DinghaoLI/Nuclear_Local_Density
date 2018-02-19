import scipy.io as sio  
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

# X = np.loadtxt("rVals_psi.txt")
# Y = np.loadtxt("zVals_psi.txt")
# Z = np.loadtxt("Psi.txt")

fig = plt.figure()

subfig1 = fig.add_subplot(2,2,1)
X = np.loadtxt("rVals_psi000.txt")
Y = np.loadtxt("zVals_psi000.txt")
Z = np.loadtxt("Psi000.txt")
surf1 = plt.contourf(X, Y, Z, 1000)
fig.colorbar(surf1)
plt.title('000')

plt.ylabel('Z')

subfig1 = fig.add_subplot(2,2,2)
X = np.loadtxt("rVals_psi001.txt")
Y = np.loadtxt("zVals_psi001.txt")
Z = np.loadtxt("Psi001.txt")
surf1 = plt.contourf(X, Y, Z, 1000)
fig.colorbar(surf1)
plt.title('001')

plt.ylabel('Z')

subfig1 = fig.add_subplot(2,2,3)
X = np.loadtxt("rVals_psi011.txt")
Y = np.loadtxt("zVals_psi011.txt")
Z = np.loadtxt("Psi011.txt")
surf1 = plt.contourf(X, Y, Z, 1000)
fig.colorbar(surf1)
plt.title('011')
plt.xlabel('R')
plt.ylabel('Z')

subfig1 = fig.add_subplot(2,2,4)
X = np.loadtxt("rVals_psi101.txt")
Y = np.loadtxt("zVals_psi101.txt")
Z = np.loadtxt("Psi101.txt")
surf1 = plt.contourf(X, Y, Z, 1000)
fig.colorbar(surf1)
plt.title('101')
plt.xlabel('R')
plt.ylabel('Z')



plt.show()