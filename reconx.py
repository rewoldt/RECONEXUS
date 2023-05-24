#!/usr/bin/env python

'''
Some simple tools to handle RECONX output.
'''
import re
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
    '''Read a null file and return a dictionary-like data object.'''

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
