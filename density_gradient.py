#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt

x, y = -np.log10(np.random.random_sample((2,10**6)))

#histogram definition
bins = [1000, 1000] # number of bins

# histogram the data
hh, locx, locy = np.histogram2d(x, y, bins=bins)

# Sort the points by density, so that the densest points are plotted last
z = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(x,y)])
idx = z.argsort()
x2, y2, z2 = x[idx], y[idx], z[idx]

plt.figure(1,figsize=(8,8)).clf()
s = plt.scatter(x2, y2, c=z2, cmap='jet', marker='.')  
plt.show()
