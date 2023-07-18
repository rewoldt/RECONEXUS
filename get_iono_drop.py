#!/usr/bin/env python3

'''
Let's grab the potential drop between two points in the ionosphere!
'''

import numpy as np
from spacepy.pybats import rim


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


def get_iono_drop(iono, x1, x2, plot=False):
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
    cpcp = iono['n_phi'].max() - iono['n_phi'].min()
    pot1 = interp(colat1, 180/np.pi*lon1)[0, 0]
    pot2 = interp(colat2, 180/np.pi*lon2)[0, 0]
    drop = max(pot1, pot2) - min(pot1, pot2)

    # Create a dictionary of returned results:
    results = {'cpcp': cpcp, 'potdrop': drop, 'pot1': pot1, 'pot2': pot2}

    if plot:
        # Use spacepy to create a quick ionosphere plot:
        fig, ax, cont, cbar = iono.add_cont('n_phi', add_cbar=True)

        # Add x1, x2 locations. Note the rotation to offset axes conventions
        # of matplotlib in polar coords:
        ax.plot([lon1+np.pi/2, lon2+np.pi/2], [colat1, colat2], 'Pw',
                mec='k', ms=12)
        ax.plot([lon1+np.pi/2, lon2+np.pi/2], [colat1, colat2], 'k--')

        # Add potential information to figure:
        fig.text(.05, .05,
                 f'CPCP={cpcp:.2f}$kV$\nNull-Null Drop={drop:.2f}$kV$')
        results['fig'] = fig
        results['ax'] = ax

        fig.tight_layout()

    return results


# Main script:
if __name__ == '__main__':
    # Open a file:
    filename = 'IE/it980504_160000_000.idl.gz'
    iono = rim.Iono(filename)

    # Hand-set our field footpoints:
    #x1 = [0.9215783971703556, -0.8818911455819862, 2.1130912688120045]
    #x2 = [0.41177270042024877, 1.1483176653635612, 2.130542854672557]
    x1 = [ 0.64182388, -0.80223776, -2.25742702]
    x2 = [ 0.63579854,  0.82509339, -2.21389861]
    x1, x2 = np.array(x1), np.array(x2)

    pots = get_iono_drop(iono, x1, x2, plot=True)
