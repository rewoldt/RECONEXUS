#!/usr/bin/env python

'''
Some tools to locate the footprints of field line traces

'''

import numpy as np

def sort_trace(X,Y,Z):
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

