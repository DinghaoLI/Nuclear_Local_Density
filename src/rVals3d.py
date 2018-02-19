import numpy as np

x,y = np.mgrid[-10:10:32j, -10:10:32j]
R = np.sqrt(x*x + y*y)

np.savetxt("rVals3d.txt", R)

print("rVals3d is saved in rVals.txt")