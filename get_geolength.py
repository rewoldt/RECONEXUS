#!/usr/bin/env python

'''
Some simple tools to fit ray traces and calculate geoeffective length
'''
import numpy as np
import matplotlib.pyplot as plt
import reconx
from scipy.stats import linregress

ray_neg_3 = reconx.read_nulls('/Users/rewoldt/RECONX/run/null_line_neg_n03_t00160000.dat')
ray_pls_3 = reconx.read_nulls('/Users/rewoldt/RECONX/run/null_line_pls_n03_t00160000.dat')

def fit_ray(R,Z):
     fit = linregress(R, Z)
     i=0
     while fit.rvalue < 0.999:
          fit = linregress(R[i:], Z[i:])
          i += 1
     return fit
     
Xn = ray_neg_3['X']
Yn = ray_neg_3['Y']
Rn = np.sqrt(Xn**2+Yn**2)
Zn = ray_neg_3['Z']

Xp = ray_pls_3['X']
Yp = ray_pls_3['Y']
Rp = np.sqrt(Xp**2+Yp**2)
Zp = ray_pls_3['Z']

fit_neg = fit_ray(Rn,Zn)
fit_pls = fit_ray(Rp,Zp)

plt.plot(Rn,Zn,Rp,Zp)
plt.plot(Rn, fit_neg.intercept + Rn*fit_neg.slope)
plt.plot(Rn, fit_pls.intercept + Rn*fit_pls.slope)
plt.show()
