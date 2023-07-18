#!/usr/bin/env python

'''
Some tools to locate the footprints of field line traces

'''

import numpy as np
import reconx

def sort_trace_old(X,Y,Z):

    ''' sort trace data if not sorted by distance '''

    data = np.concatenate((X[:, np.newaxis],
                       Y[:, np.newaxis],
                       Z[:, np.newaxis]),
                       axis=1)
    distances = np.linalg.norm(data, axis=1)
    # Get sorted indices
    sorted_indices = np.argsort(distances)
    # Sort the data based on the indices
    sorted_data = data[sorted_indices]
    return sorted_data

def sort_trace(null_file):

    ''' sort trace data if not sorted by distance '''

    X = null_file['X']
    Y = null_file['Y']
    Z = null_file['Z']

    data = np.concatenate((X[:, np.newaxis],
                       Y[:, np.newaxis],
                       Z[:, np.newaxis]),
                       axis=1)
    distances = np.linalg.norm(data, axis=1)
    # Get sorted indices
    sorted_indices = np.argsort(distances)
    # Sort the data based on the indices
    sorted_data = data[sorted_indices]

    X_new, Y_new, Z_new = sorted_data[:,0],sorted_data[:,1],sorted_data[:,2]
    print (X_new, Y_new, Z_new)

    return X_new, Y_new, Z_new

def get_footprints(sorted_data):

    ''' get the footprints corresponding to each null '''

    return sorted_data[0].tolist()

if __name__ == '__main__':

    ray_neg = reconx.read_nulls('/Users/rewoldt/RECONX/run/null_line_neg_s03_111_t00160000.dat')
    ray_pls = reconx.read_nulls('/Users/rewoldt/RECONX/run/null_line_pls_s03_125_t00160000.dat')

    Xn = ray_neg['X']
    Yn = ray_neg['Y']
    Zn = ray_neg['Z']
    Xp = ray_pls['X']
    Yp = ray_pls['Y']
    Zp = ray_pls['Z']
    
    Xp_sorted,Yp_sorted,Zp_sorted = sort_trace(ray_pls)
    Xn_sorted,Yn_sorted,Zn_sorted = sort_trace(ray_neg)
    #Pp = sort_trace(Xp,Yp,Zp)
    #Pn = sort_trace(Xn,Yn,Zn)
    #coord_pls = get_footprints(Pp)
    #coord_neg = get_footprints(Pn)
    print ('Last entry for negative trace:', Xn[-1],Yn[-1],Zn[-1], '\n', 'First entry of sorted negative null trace:', Xn_sorted[0],Yn_sorted[0],Zn_sorted[0])
    print ('Last entry for positive trace:',Xp[-1],Yp[-1],Zp[-1], '\n', 'First entry of sorted negative null trace:', Xp_sorted[0],Yp_sorted[0],Zp_sorted[0])