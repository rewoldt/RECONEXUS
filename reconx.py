#!/usr/bin/env python

'''
Some simple tools to handle RECONX output.
'''
import re
import os
from glob import glob
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt

from spacepy.pybats import ImfInput, bats, rim
# Important constants:
RE = 6371000  # in meters!


def dot_product(u, v):
    '''
    Calculate and return $u \cdot v$ for inputs u, v.

    Paramters
    =========
    u, v : 3-element numpy vectors
        Vectors to use for cross product.

    Returns
    =======
        Dot product result.
    '''
    return sum([un*vn for un, vn in zip(u, v)])


def cross_product(u, v):
    '''
    Calculate and return $u \times v$ for inputs u, v.

    Parameters
    =========
    u, v : 3-element numpy vectors
        Vectors to use for cross product.

    Returns
    =======
    w : 3-element numpy vector
        Cross product result.
    '''

    w = np.array([u[1]*v[2] - u[2]*v[1],
                  u[2]*v[0] - u[0]*v[2],
                  u[0]*v[1] - u[1]*v[0]])

    return w

def map_from_IB(lat, R=2.5):
    '''
    Given a geomagnetic latitude, `lat` (radians) at some distance `R`
    at the MHD inner boundary, map down to ionospheric altitudes (R~=1)
    assuming a dipole configuration.
    '''

    return np.arccos(np.sqrt(np.cos(lat)**2/R))


def xyz_to_lonlat(x, domap=True):
    '''
    Given XYZ as a 3-element numpy array, convert to lon and lat; return
    lon and lat in radians.

    If `domap` is True, the latitude will be mapped from its current radius
    down to R=1 assuming magnetic dipole geometry.
    '''

    R = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    xy = np.sqrt(x[0]**2 + x[1]**2)

    lon, lat = np.arctan2(x[1], x[0]), np.arctan2(x[2], xy)

    if R > 1 and domap:
        lat = map_from_IB(lat, R=R)

    return lon, lat


def get_iono_drop(iono, x1, x2):
    '''
    Given a `pybats.rim.Iono` object, `iono`, and two points (either in
    SM XYZ or SM lon/lat), determine the electric potential drop between
    those points.

    `x1` and `x2` should be either 3-element numpy arrays (for SM cartesian
    XYZ in units of RE) or 2-element numpy arrays (for SM lon/lat in radians).

    If `plot` is set to True, a figure is created showing the locations of
    `x1` and `x2` in the northern ionosphere.

    A dictionary of results is returned to users which contains the potential
    values at x1, x2; the potential drop between them; the total CPCP; and
    any plot objects if `plot` is set to True.
    '''

    from scipy.interpolate import RectBivariateSpline as Spline

    # If in cartesian XYZ, convert to SM lat/lon:
    if x1.shape != 3:
        lon1, lat1 = xyz_to_lonlat(x1)
    else:
        lon1, lat1 = x1[:]

    if x2.shape != 3:
        lon2, lat2 = xyz_to_lonlat(x2)
    else:
        lon2, lat2 = x2[:]

    # Positive longitudes only:
    lon1 += 2*np.pi*(lon1 < 0)
    lon2 += 2*np.pi*(lon2 < 0)

    # Lat to colat in degrees:
    colat1 = 90 - 180/np.pi * lat1
    colat2 = 90 - 180/np.pi * lat2

    # Create interpolator object:
    lons, lats = iono['n_psi'][0, :], iono['n_theta'][:, 0]
    interp = Spline(lats, lons, iono['n_phi'])

    # Get potentials:
    pot1 = interp(colat1, 180/np.pi*lon1)[0, 0]
    pot2 = interp(colat2, 180/np.pi*lon2)[0, 0]
    drop = max(pot1, pot2) - min(pot1, pot2)

    return drop

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


def read_nulls(filename, reorder=False):
    '''
    Read a null file and return a dictionary-like data object.

    Parameters
    ----------
    filename : str
        Name of null file to open.
    reorder : bool, defaults to False
        Reorder the  points in the file such that they represent, in order,
        a trace along a magnetic field line.

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

    # Re-arrange points:
    if reorder:
        # Find the spot where the trace flips direction:
        r = np.sqrt(data['X']**2 + data['Y']**2 + data['Z']**2)
        i = np.argmax(np.abs(r[1:] - r[:-1]))

        # Re-sort the lines:
        data['X'] = np.append(data['X'][-1:i+1:-1], data['X'][:i])
        data['Y'] = np.append(data['Y'][-1:i+1:-1], data['Y'][:i])
        data['Z'] = np.append(data['Z'][-1:i+1:-1], data['Z'][:i])

    return data


class NullPair(dict):
    '''
    Class for handling a single Null Pair.

    Parameters
    ----------
    x_pos, x_neg : dmarray
        The XYZ position of the positive and negative nulls, respectively
    inull : int
        The number of the null pair.
    usw, bsw : 3-element numpy vectors
        The solar wind velocity and magnetic field corresponding to the time
        that the null was found. Units should be km/s and nT, respectively.
    cpcp : int, defaults to None
        Cross polar cap potential (cpcp) corresponding to point in time when
        null found.
    potdrop : int, defaults to None
        Potential drop between fotpoints in the ionosphere corresponding to
        point in time when null found.
    path : str, defaults to empty string
        Path to file location, used for opening line files.
    time : datetime.datetime, defaults to None
        Date time corresponding to point in time when null found.
    iline : int, defaults to 1
        If there are more than 1 lines associated with this null, this
        kwarg sets which to read.
    '''

    def __init__(self, x_pos, x_neg, inull, usw, bsw, ionofile='', path='',
                 time=None, iline=1):
        # Initialize as a dictionary:
        super(NullPair, self).__init__()

        # Set positions:
        self['x_pos'], self['x_neg'] = x_pos, x_neg
        self['time'] = time
        self['inull'] = inull

        # Save solar wind conditions:
        self.u, self.b = usw, bsw

        # Get magnetic field line info for each null.
        f_pos = path + f'null_line_pls_n{inull:02d}_{iline:03d}' + \
            f'_e{time:%Y%m%d-%H%M%S}.dat'
        f_neg = path + f'null_line_neg_n{inull:02d}_{iline:03d}' + \
            f'_e{time:%Y%m%d-%H%M%S}.dat'
        self['linefiles'] = [f_pos, f_neg]

        # Are there line files for this null? If not, crash.
        if not os.path.exists(f_pos) or not os.path.exists(f_neg):
            # Line tracing failed for this null. Fill defaults and stop.
            self['posline'] = False
            self['negline'] = False
            return

        # Open field line file and read.
        self['posline'] = read_nulls(f_pos, reorder=True)
        self['negline'] = read_nulls(f_neg, reorder=True)

        self.calc_geopot()
        self.get_cpcp(ionofile)
        self.calc_potdrop(ionofile)


    def calc_geopot(self):
        '''
        Calculate potential drop across geoeffective length, save as
        self['geopot']. Units are kV.
        '''

        # Get start and end point of integration:
        s1 = np.array([self['posline']['X'][10],
                       self['posline']['Y'][10],
                       self['posline']['Z'][10]])
        s2 = np.array([self['negline']['X'][10],
                       self['negline']['Y'][10],
                       self['negline']['Z'][10]])
        s = RE * (s2 - s1)  # Integration path as a vector in SI units (m)]

        # Perform integration assuming constant values across line.
        # Unit conversion is nT->T; result is in kV.
        self['geopot'] = 1E-9*dot_product(cross_product(self.u, self.b), s)


    def get_cpcp(self, ionofile):
        '''
        Open `ionofile` and calculate cpcp associated with null pair.
        '''
        # Can't find file? No filename given? Set defaults and bail.

        if os.path.exists(ionofile):
            self.iono = rim.Iono(ionofile)
            self['cpcp'] = self.iono['n_phi'].max() - self.iono['n_phi'].min()
        else: self['cpcp'] = 0
        return

    def calc_potdrop(self, ionofile):
        '''
        Calculate potential drop across footpoints in the ionosphere,
        save as self['potdrop']. Units are kV.

        '''
        # Get start and end point of integration:
        x1 = np.array([self['posline']['X'][-1],
                       self['posline']['Y'][-1],
                       self['posline']['Z'][-1]])
        x2 = np.array([self['negline']['X'][-1],
                       self['negline']['Y'][-1],
                       self['negline']['Z'][-1]])
        if os.path.exists(ionofile):
            self.iono = rim.Iono(ionofile)
            self['potdrop'] = get_iono_drop(self.iono, x1, x2)
        else: self['potdrop'] = 0
        return

    def __repr__(self):
        return f'NullPair at [{self["x_pos"]}, {self["x_neg"]}]'


class NullGroup(list):
    '''
    Class for handling magnetic null pairs and all associated information.

    To instantiate,

    Paramters
    =========
    nullfile : str
        Path to null file, e.g., "NegNulls_e20150321-040000.dat". These files
        contain a list of nulls taken at a certain time.
    rundir : str, defaults to None
        If given, associates a SWMF run directory with this null set.
    imffile : str, defaults to None
        If given, associates an SWMF IMF input file with this null set.
        Will be used to set solar wind conditions when calculating potential
        drops across lines.
    '''

    def __init__(self, nullfile, rundir='', imffile=None):
        # Initialize as a dictionary:
        super(NullGroup, self).__init__()

        # Get first null file, find match.
        if 'Plus' in nullfile:
            f_plus = nullfile
            f_nega = nullfile.replace('Plus', 'Neg')
        elif 'Neg' in nullfile:
            f_nega = nullfile
            f_plus = nullfile.replace('Neg', 'Plus')

        # Get the path of the files:
        path = '/'.join(nullfile.split('/')[:-1])
        if path:
            path = path + '/'

        # Open null files and start mining information.
        posn, negn = read_nulls(f_plus), read_nulls(f_nega)

        # Get datetime from file names.
        self.time = dt.datetime.strptime(nullfile[-19:-4], '%Y%m%d-%H%M%S')

        # Find IMF file, set U and B in solar wind.
        # If nothing found, set default U, B conditions.
        if not imffile:  # No imffile set?
            # Try to auto-find:
            #files = glob(rundir + 'IMF*.dat')
            files = glob(rundir + 'imf*.dat')
            if files:
                imffile = files[0]
        # Set IMF conditions (use default of no imffile found.)
        self.get_imf(imffile)

        files = glob(rundir+f'IE/it{nullfile[-17:-11]}_{nullfile[-10:-4]}_000.*')
        if files:
            ionofile = files[0]
        else: ionofile=''

        #ionofile = glob(rundir+f'IE/it{self.time: %Y%m%d-%H%M%S}_000.idl')[0]
        # Set cpcp from ionosphere files.
        #self.get_cpcp(ionofile)

        # For each null, create a NullPair object and store.
        self.nnulls = min(posn['X'].size, negn['X'].size)
        #self.nnulls = posn['X'].size
        for i in range(self.nnulls):
            x_pos = np.array([posn['X'][i], posn['Y'][i], posn['Z'][i]])
            x_neg = np.array([negn['X'][i], negn['Y'][i], negn['Z'][i]])
            self.append(NullPair(x_pos, x_neg, i+1, self.u, self.b, ionofile=ionofile, time=self.time, path=path))

    def get_imf(self, imffile, defaultu=np.array([-400, 0, 0]),
                defaultb=np.array([0, 0, -10])):
        '''
        Open `imffile` and interpolate to self.time. If imffile does not
        exist, then default to `defaultu`, `defaultb`.

        If `imffile` is a spacepy.pybats.ImfInput object, use that.
        If `imffile` is a string, open the file as a ImfInput object.
        '''

        from matplotlib.dates import date2num

        if type(imffile) is str:
            # Can't find file? No filename given? Set defaults and bail.
            if not imffile or not os.path.exists(imffile):
                self.imf = None
                self.u = defaultu
                self.b = defaultb
                return
            # Open the file if found:
            self.imf = ImfInput(imffile)

        elif type(imffile) is ImfInput:
            self.imf = imffile
        else:
            raise ValueError('Unrecognized type for imffile')

        t_imf = date2num(self.imf['time'])
        t_now = date2num(self.time)

        self.u, self.b = np.zeros(3), np.zeros(3)
        self.u[0] = np.interp(t_now, t_imf, self.imf['ux'])
        self.u[1] = np.interp(t_now, t_imf, self.imf['uy'])
        self.u[2] = np.interp(t_now, t_imf, self.imf['uz'])
        self.b[0] = np.interp(t_now, t_imf, self.imf['bx'])
        self.b[1] = np.interp(t_now, t_imf, self.imf['by'])
        self.b[2] = np.interp(t_now, t_imf, self.imf['bz'])

        return