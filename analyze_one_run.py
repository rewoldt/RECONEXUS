#!/usr/bin/env python3

'''
Tools required to analyze a single simulation.

Reqs:
    One SWMF run directory (variable `rundir_test`)
    One RECONEXUS directory (variable `directory_test`)
    The associated IMF file (variable `imffile`)

Process:
    Set the variables above.
    Run "process_simulation" to create a timeseries files containing all
    values.
    Plot to your heart's content.

Underlying Assumptions:
    Reconexus was run using the separator method; only endpoints are
    considered.
    Only set of trace files per timestep exists.
    The ionospheric pedersen conductance has been manually entered into the
    IMF input file.
'''

import os
import datetime as dt
from glob import glob
import pickle

import numpy as np
import matplotlib.pyplot as plt
import spacepy.plot as splot
from spacepy.pybats import ImfInput

import reconx


# Use Spacepy's plot style:
splot.style()

# Simulation locations:
runs = {}

# Add runs to this prefix:
prefix = "/mnt/y/rewoldt/rxn_runs/"
runs['vsw'] = {
    'directory_test': prefix + "reconx_data/run_DeltaVsw_surface_long/",
    'rundir_test': prefix + "rxn_DeltaVsw_new/",
    'imffile': prefix + "rxn_DeltaVsw_new/imf_mf_DeltaVsw_by.dat"
}

runs['bz'] = {
    'directory_test': prefix + "reconx_data/run_deltaBz_long_surface/",
    'rundir_test': prefix + "rxn_deltaB_long/",
    'imffile': prefix + "rxn_deltaB_long/imf_mf_DeltaBz_by.dat"
}

# NUMERICAL CONSTANTS:
ccm2SI = 1.6726E-27*1.0E6  # Conversion from #/cm-3 to kg/m-3
mu0 = 4*np.pi*1.0E-7  # Permeability of free space


def calc_DeltaV(geo_pot, bsw, rho, Sigma_P, debug=False):
    '''
    Calculate theoretical CPCP from Kivelson and Ridley.
    '''

    # Unit conversions:
    rhosw = rho * ccm2SI  # conversion for
    bsw = bsw*1.0E-9

    vA = bsw/np.sqrt(mu0*rhosw)  # V_A = B / (mu*rho)^(1/2) m/s
    Sigma_A = 1/(mu0*vA)

    DeltaV = (abs(geo_pot)*Sigma_A)/(Sigma_P+Sigma_A)
    if debug:
        print(f'alfven conductance={Sigma_A};' +
              f'Pederson Conductance {Sigma_P};' +
              f'alfven speed={vA}; DeltaV={DeltaV}')

    return Sigma_A, DeltaV


def process_simulation(runname):
    '''
    Calculate various potentials, etc and place into a single timeseries file.

    Resulting data will be a dictionary with time, runtime (hours), and
    associated values.
    '''

    # Load IMF file:
    imf = ImfInput(runs[runname]['imffile'])

    # Set corresponding directories;
    dir_test = runs[runname]['directory_test']
    dir_run = runs[runname]['rundir_test']

    # Get list of null point files:
    files = glob(dir_test + 'NegNulls*.dat')
    files.sort()  # Keep things in order!

    # Create a nice output data object:
    data = {}
    outvars = ['time', 'cpcp', 'geopot', 'ifpot', 'gel', 'alf',
               'krpot', 'bz', 'b', 'u', 'ux', 'n', 'ped', 'sigmaA']
    for v in outvars:
        data[v] = []

    for f in files:
        # Try to open the NullGroup; it may fail if there is no matching
        # trace file. Skip if this is true.
        null = reconx.NullGroup(f, dir_run, imf)

        # Make sure we just have 1 null pair in our group.
        if len(null) > 1:
            print(f'WARNING: MULTIPLE NULLS AT {null.time}')
            continue

        if len(null) < 1:
            print(f'WARNING: ZERO NULLS AT {null.time}')
            continue

        if not null[0]['posline']:
            print(f"Null TRACING not found for {f}. Skipping...")
            continue

        # Save all info into our data array.
        data['time'].append(null.time)
        data['b'].append(null.btot)
        data['bz'].append(null.b[-1])
        data['u'].append(null.utot,)
        data['ux'].append(null.u[0])
        data['n'].append(null.n)
        data['alf'].append(null.alf)
        data['ped'].append(null.ped)

        data['geopot'].append(null[0]["geopot"])
        data['cpcp'].append(null[0]["cpcp"])
        data['ifpot'].append(null[0]["potdrop"])
        data['gel'].append(null[0]["gel"])
        sigA, krpot = calc_DeltaV(null[0]["geopot"], null.btot,
                                  null.n, null.ped)
        data['krpot'].append(krpot)
        data['sigmaA'].append(sigA)

    # Convert lists to numpy arrays:
    for v in outvars:
        data[v] = np.array(data[v])
    # Calculate run time in hours:
    tstart = data['time'][0]
    data['runtime'] = np.array([(t - tstart).total_seconds()/3600.
                                for t in data['time']])

    # Get absval of Geopot, which is sometimes negative.
    data['geopot'] = np.abs(data['geopot'])

    with open(f'reconexus_results_{runname}.pkl', 'wb') as f:
        pickle.dump(data, f)

    return data


def plot_summary(runname):
    '''
    Given a results pickle, load the data and create a quick-look summary
    plot.
    '''

    # Load a pickle.
    with open(f'reconexus_results_{runname}.pkl', 'rb') as f:
        data = pickle.load(f)

    fig, (a1, a2, a3, a4, a5) = plt.subplots(5, 1, sharex=True, figsize=(9, 9))
    fig.suptitle(f'Results from {runname.capitalize()} Investigation')

    a1.plot(data['runtime'], data['bz'], 'o-')
    a1twin = a1.twinx()
    a1twin.plot(data['runtime'], data['n'], 'ro-')
    a1.set_ylabel('IMF Bz', c='C0')
    a1twin.set_ylabel(r'$\rho_{sw}$', c='r')

    a2.plot(data['runtime'], data['u'], 'o-')
    a2.set_ylabel(r'$V_{SW}$', c='C0')
    a2twin = a2.twinx()
    a2twin.plot(data['runtime'], data['alf'], 'o-', c='C1')
    a2twin.set_ylabel(r'$V_{Alf, SW}$', c='C1')

    a3.plot(data['runtime'], data['ped'], label=r'$\Sigma_{Ped,Iono}$')
    a3.plot(data['runtime'], data['sigmaA'], label=r'$\Sigma_{SW}$')
    a3.legend(loc='best')
    a3.set_ylabel('Siemens')

    a4.plot(data['runtime'], data['gel'], 'o-')
    a4.set_ylabel('GEL ($km$)')

    a5.plot(data['runtime'], data['geopot'], 'o-', label=r'Reconexus $\Phi$')
    a5.plot(data['runtime'], data['cpcp'], 'o-', label='CPCP')
    a5.plot(data['runtime'], data['krpot'], 'o-', label='K&R Pot')
    a5.plot(data['runtime'], data['ifpot'], 'o-', label='Footpoint Pot')
    a5.legend(loc='best')
    a5.set_ylabel(r'$\Phi$ ($kV$)')

    fig.tight_layout()


def compare_gel_2times(runname, t1=dt.datetime(1998, 5, 4, 8, 0, 0),
                       t2=dt.datetime(1998, 5, 5, 11, 0, 0)):
    '''
    Examine GEL dynamics between two times.

    t1 and t2 should be datetimes. :P
    '''

    # Open the files we want:
    # Load a pickle with reconnexus stuff:
    with open(f'reconexus_results_{runname}.pkl', 'rb') as f:
        data = pickle.load(f)

    # Shortcut vars:
    dirrxn = runs[runname]['directory_test']

    # Open our RxN lines:
    strtime = f"{t1:%Y%m%d-%H%M%S}"
    lp1 = reconx.read_nulls(dirrxn + f"null_line_pls_n01_001_e{strtime}.dat")
    ln1 = reconx.read_nulls(dirrxn + f"null_line_neg_n01_001_e{strtime}.dat")

    strtime = f"{t2:%Y%m%d-%H%M%S}"
    lp2 = reconx.read_nulls(dirrxn + f"null_line_pls_n01_001_e{strtime}.dat")
    ln2 = reconx.read_nulls(dirrxn + f"null_line_neg_n01_001_e{strtime}.dat")


    fig = plt.figure(figsize=[16, 12])
    a1, a2 = fig.add_subplot(2, 2, 1), fig.add_subplot(2, 2, 2, projection='3d')
    a3 = fig.add_subplot(2, 1, 2)

    a1.plot(lp1['X'], lp1['Y'], 'ro', ln1['X'], ln1['Y'], 'ro', label='Time 1')
    a1.plot(lp2['X'], lp2['Y'], 'bo', ln2['X'], ln2['Y'], 'bo', label='Time 2')
    a1.legend(loc='best')

    a2.plot(lp1['X'], lp1['Y'], lp1['Z'], '.r')
    a2.plot(ln1['X'], ln1['Y'], ln1['Z'], '.r', label='Time 1')
    a2.plot(lp2['X'], lp2['Y'], lp2['Z'], '.b')
    a2.plot(ln2['X'], ln2['Y'], ln2['Z'], '.b', label='Time 2')

    a3.plot(data['runtime'], data['geopot'], 'o-', label=r'Reconexus $\Phi$')
    a3.plot(data['runtime'], data['cpcp'], 'o-', label='CPCP')
    a3.plot(data['runtime'], data['krpot'], 'o-', label='K&R Pot')
    a3.plot(data['runtime'], data['ifpot'], 'o-', label='Footpoint Pot')
    a3.legend(loc='best')
    a3.set_ylabel(r'$\Phi$ ($kV$)')
