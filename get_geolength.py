#!/usr/bin/env python                                                                                                       

'''                                                                                                                          
Some simple tools to fit ray traces and calculate geoeffective length                                                        
'''

import numpy as np
import matplotlib.pyplot as plt
import reconx
import get_footprints

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

def get_potential(null_file_neg, null_file_pls,n, plot=False):

     Xp_sorted,Yp_sorted,Zp_sorted = get_footprints.sort_trace(null_file_pls)
     Xn_sorted,Yn_sorted,Zn_sorted = get_footprints.sort_trace(null_file_neg)

     Xn = Xn_sorted[n:]
     Yn = Yn_sorted[n:]
     Zn = Zn_sorted[n:]
     Xp = Xp_sorted[n:]
     Yp = Yp_sorted[n:]
     Zp = Zp_sorted[n:]

     RE2km = 6378.14 *(10**3) #m

     #geo_potential = []
     dist = []
     length = []
     length_dict = {}

     uxb = calc_uxb(u,B)

     for i in range(len(Xn)):
          for j in range(len(Xp)):
               distance = (np.sqrt((Xn[i]-Xp[j])**2 + (Yn[i]-Yp[j])**2 + (Zn[i]-Zp[j])**2)*RE2km)
               dist.append(distance)
               L = (Xn[i]-Xp[j])*RE2km, (Yn[i]-Yp[j])*RE2km, (Zn[i]-Zp[j])*RE2km
               length.append([(Xn[i]-Xp[j])*RE2km, (Yn[i]-Yp[j])*RE2km, (Zn[i]-Zp[j])*RE2km])
               #geo_potential.append(dot_product(uxb,L))

     min_distance_index = np.argmin(dist)
     min_distance = dist[min_distance_index]

     print (min_distance_index, min_distance, length[min_distance_index], length[min_distance_index][0], length[min_distance_index][1], length[min_distance_index][2], (np.sqrt((length[min_distance_index][0])**2 + (length[min_distance_index][1])**2 + (length[min_distance_index][2])**2)))

     geo_potential = dot_product(uxb,length[min_distance_index])

     if plot:
               
          fig = plt.figure()
          ax = fig.add_subplot(projection='3d')

          ax.scatter3D(0,0,0, s=50,color='Black')
          ax.scatter3D(Xn*RE2km, Yn*RE2km, Zn*RE2km, color = 'Blue', label='Negative trace')
          ax.scatter3D(Xp*RE2km, Yp*RE2km, Zp*RE2km, color = 'Red', label='Positive trace')
          ax.set_xlabel('X')
          ax.set_ylabel('Y')
          ax.set_zlabel('Z')
          ax.quiver(Xp*RE2km, Yp*RE2km, Zp*RE2km,length[min_distance_index][0], length[min_distance_index][1], length[min_distance_index][2], label='uxB vector')
          #for i in range(len(Xn)):
          #     for j in range(len(Xp)):
          #          ax.plot([Xn[i],Xp[j]], [Yn[i], Yp[j]], [Zn[i], Zp[j]])
          ax.legend()
          #ax.scatter3D(length)
          plt.show()

     return dist, geo_potential

if __name__ == '__main__':
     ray_neg = reconx.read_nulls('/Users/rewoldt/RECONX/run/null_line_neg_s03_111_t00160000.dat')
     ray_pls = reconx.read_nulls('/Users/rewoldt/RECONX/run/null_line_pls_s03_125_t00160000.dat')

     ux = -450*10**3 #m/s
     uy = 0
     uz = 0

     Bx = 0
     By = 2*10**(-9) #T
     Bz = -10*10**(-9) #T
     #Bz = 5*(10**(-9)) #T

     B = [Bx, By, Bz]
     u = [ux,uy,uz]

     dist, geo_potential = get_potential(ray_neg, ray_pls,-30,plot=True)

     print ('geo length: ', min(dist), 'm', 'geo potential', geo_potential, 'V', 'Should be getting ~95.72 kV')