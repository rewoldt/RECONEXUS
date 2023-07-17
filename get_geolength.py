#!/usr/bin/env python                                                                                                       

'''                                                                                                                          
Some simple tools to fit ray traces and calculate geoeffective length                                                        
'''

import numpy as np
import matplotlib.pyplot as plt
import reconx

def cross_product(u,v):  
    dim = len(u)
    s = []
    for i in range(dim):
        if i == 0:
            j,k = 1,2
            s.append(u[j]*v[k] - u[k]*v[j])
        elif i == 1:
            j,k = 2,0
            s.append(u[j]*v[k] - u[k]*v[j])
        else:
            j,k = 0,1
            s.append(u[j]*v[k] - u[k]*v[j])
    return s

def dot_product(u,v):
    return sum([un*vn for un,vn in zip(u,v)])

def calc_uxb(u,B):
    uxb = cross_product(u,B)
    return uxb

def get_potential(null_file_neg, null_file_pos):

    Xn = null_file_neg['X'][:200]
    Yn = null_file_neg['Y'][:200]
    Zn = null_file_neg['Z'][:200]
    Xp = null_file_pos['X'][:200]
    Yp = null_file_pos['Y'][:200]
    Zp = null_file_pos['Z'][:200]

    RE2km = 6378.14*(10**3) #m

    geo_potential = []
    dist = []
    length = []

    uxb = calc_uxb(u,B)

    for i in range(len(Xn)):
        for j in range(len(Xp)):
            dist.append(np.sqrt((Xn[i]-Xp[j])**2 + (Yn[i]-Yp[j])**2 + (Zn[i]-Zp[j])**2)*RE2km)
            L = (Xn[i]-Xp[j])*RE2km, (Yn[i]-Yp[j])*RE2km, (Zn[i]-Zp[j])*RE2km
            length.append([(Xn[i]-Xp[j])*RE2km, (Yn[i]-Yp[j])*RE2km, (Zn[i]-Zp[j])*RE2km])
            print (L, uxb)
            print (np.dot(uxb,L), dot_product(uxb,L))
            geo_potential.append(dot_product(uxb,L))

    return dist, min(geo_potential)


if __name__ == '__main__':
    ray_neg = reconx.read_nulls('/Users/rewoldt/RECONX/run/null_line_neg_s03_111_t00160000.dat')
    ray_pls = reconx.read_nulls('/Users/rewoldt/RECONX/run/null_line_pls_s03_125_t00160000.dat')

    ux = -450*(10**3) #km/s
    uy = 0
    uz = 0

    Bx = 0
    By = 2*(10**(-9)) #nT
    Bz = -10*(10**(-9)) #nT

    B = [Bx, By, Bz]
    u = [ux,uy,uz]

    dist, geo_potential = get_potential(ray_neg, ray_pls)

    print ('geo length: ', min(dist), 'geo potential', geo_potential, 'Should be getting ~95.72')

