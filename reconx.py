#!/usr/bin/env python

'''
Some simple tools to handle RECONX output.
'''
import re
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt


def read_separator(filename):
    '''Read a separator file and return a dictionary-like data object.'''

    with open(filename, 'r') as f:
        # Read one line, parse header:
        line = f.readline()
        head = re.findall('\"(.+?)\s\[(.+?)\]\"', line)

        # Extract variable names and units.
        var, unit = [], []
        for pair in head:
            var.append(pair[0])
            unit.append(pair[1])

        # Skip next line in header:
        f.readline()

        # Read remainder of lines:
        lines = f.readlines()

    # Create container for data:
    data = {}
    for v in var:
        data[v] = np.zeros(len(lines))

    # Put data into the container
    for i, l in enumerate(lines):
        parts = l.split()
        for v, p in zip(var, parts):
            data[v][i] = p

    return data


def read_nulls(filename):
    '''
    Read a null file and return a dictionary-like data object.

    Parameters
    ----------
    filename : str
        Name of null file to open.

    Returns
    -------
    data : dictionary
        Dictionary of nulls and null information.
    '''

    with open(filename, 'r') as f:
        # Read one line, parse header:
        line = f.readline()
        head = re.findall('\"(.+?)\"', line)

        # Extract variable names and units.
        var, unit = [], []
        for pair in head:
            var.append(pair[0])
            #unit.append(pair[1])

        # Skip next line in header:
        f.readline()

        # Read remainder of lines:
        lines = f.readlines()

    # Create container for data:
    data = {}
    for v in var:
        data[v] = np.zeros(len(lines))

    # Put data into the container
    for i, l in enumerate(lines):
        parts = l.split()
        for v, p in zip(var, parts):
            data[v][i] = p

    return data


def read_separator_plane(filename):
    '''Read a separator file from plane method and return a dictionary-like data object.'''

    with open(filename, 'r') as f:
        # Read one line, parse header:
        line = f.readline()
        head = re.findall('\"(.+?)\"', line)

        # Extract variable names and units.
        var, unit = [], []
        for pair in head:
            var.append(pair[0])
            #unit.append(pair[1])

        # Skip next line in header:
        f.readline()

        # Read remainder of lines:
        lines = f.readlines()

    # Create container for data:
    data = {}
    for v in var:
        data[v] = np.zeros(len(lines))

    # Put data into the container
    for i, l in enumerate(lines):
        parts = l.split()
        for v, p in zip(var, parts):
            data[v][i] = p

    return data


def read_sphere(filename):
    '''Read a separator file and return a dictionary-like data object.'''

    with open(filename, 'r') as f:
        # Read one line, parse header:
        line = f.readline()
        head = re.findall('\"(.+?)\"', line)

        # Extract variable names and units.
        var, unit = [], []
        for pair in head:
            var.append(pair[0])
        #    unit.append(pair[1])

        # Skip next line in header:
        f.readline()

        # Read remainder of lines:
        lines = f.readlines()

    # Create container for data:
    data = {}
    for v in var:
        data[v] = np.zeros(len(lines))

    # Put data into the container
    for i, l in enumerate(lines):
        parts = l.split()
        for v, p in zip(var, parts):
            data[v][i] = p

    return data


class NullPair(dict):
    '''
    Class for handling a single Null Pair.

    Parameters
    ----------
    x_pos, x_neg : dmarray
        The XYZ position of the positive and negative nulls, respectively
    '''

    def __init__(self, x_pos, x_neg, time, inull, iline=1):
        # Initialize as a dictionary:
        super(NullPair, self).__init__()

        # Set positions:
        self['x_pos'], self['x_neg'] = x_pos, x_neg

        # Get magnetic field line info for each null.
        f_pos = f'null_line_neg_o{inull:02d}_{iline:03d}' + \
                f'_e{time:%Y%m%d-%H%M%S}.dat'

        # Open field line file and read.
        b = read_nulls(f_pos)


    def __repr__(self):
        return f'NullPair at [{self["x_pos"]}, {self["x_neg"]}]'


class NullGroup(list):
    '''
    Class for handling magnetic null pairs and all associated information.

    To instantiate,
    '''

    def __init__(self, nullfile, rundir):
        # Initialize as a dictionary:
        super(NullGroup, self).__init__()

        # Get first null file, find match.
        if 'Plus' in nullfile:
            f_plus = nullfile
            f_nega = nullfile.replace('Plus', 'Neg')
        elif 'Neg' in nullfile:
            f_nega = nullfile
            f_plus = nullfile.replace('Neg', 'Plus')

        # Open null files and start mining information.
        posn, negn = read_nulls(f_plus), read_nulls(f_nega)

        # Get datetime from file names.
        self.time = dt.datetime.strptime(nullfile[-19:-4], '%Y%m%d-%H%M%S')

        # For each null, create a NullPair object and store.
        self.nnulls = posn['X'].size
        for i in range(self.nnulls):
            x_pos = np.array([posn['X'][i], posn['Y'][i], posn['Z'][i]])
            x_neg = np.array([negn['X'][i], negn['Y'][i], negn['Z'][i]])
            self.append(NullPair(x_pos, x_neg, self.time))
