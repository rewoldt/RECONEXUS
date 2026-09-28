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

import datetime as dt
from glob import glob
import pickle

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import spacepy.plot as splot
from spacepy.pybats import ImfInput, bats

import reconx


# Use Spacepy's plot style:
splot.style()

# Simulation locations:
runs = {}

# Add runs to this prefix:
prefix = "/data/rewoldt/rxn_runs/"
runs['vsw'] = {
    'rxn_dir': prefix + "reconx_data/run_DeltaVsw_surface_long/",
    'swmf_dir': prefix + "rxn_DeltaVsw_long/",
    'imffile': prefix + "rxn_DeltaVsw_long/imf_mf_DeltaVsw_by.dat"
}

runs['bz'] = {
    'rxn_dir': prefix + "reconx_data/run_DeltaBz_long_null/",
    'swmf_dir': prefix + "rxn_DeltaBz_long/",
    'imffile': prefix + "rxn_DeltaBz_long/imf_mf_DeltaBz_by.dat"
}
    
runs['n'] = {
    'rxn_dir': prefix + "reconx_data/run_DeltaN_long_null/",
    'swmf_dir': prefix + "rxn_DeltaN_long/",
    'imffile': prefix + "rxn_DeltaN_long/imf_mf_DeltaN_by.dat"
}

runs['cond'] = {
    'rxn_dir': prefix + "reconx_data/run_cond_null/",
    'swmf_dir': prefix + "rxn_conductance/",
    'imffile': prefix + "rxn_conductance/imf_mf_cond_by.dat"
                    }

runs['bz_vsw750'] = {
    'rxn_dir': prefix + "../RECONX/run_DeltaBz_Vsw750_xlim20_null/",
    'swmf_dir': prefix + "rxn_DeltaBz_Vsw750/",
    'imffile': prefix + "rxn_DeltaBz_Vsw750/imf_DeltaBz_Vsw750_by.dat"
                    }

runs['bz_20_single'] = {
    'rxn_dir': prefix + "../RECONX/run_DeltaBz20_single_null/",
    'swmf_dir': prefix + "rxn_DeltaBz20_single/",
    'imffile': prefix + "rxn_DeltaBz20_single/imf_mf_bzturn_by.dat"
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
    dir_rxn = runs[runname]['rxn_dir']
    dir_swmf = runs[runname]['swmf_dir']

    # Get list of null point files:
    files = glob(dir_rxn + 'NegNulls*.dat')
    files.sort()  # Keep things in order!

    # Create a nice output data object:
    data = {}
    outvars = ['time', 'cpcp', 'geopot', 'ifpot', 'gel', 'alf',
               'krpot', 'bz', 'b', 'u', 'ux', 'n', 'ped', 'sigmaA',
               'pram']
    for v in outvars:
        data[v] = []

    for f in files:
        # Try to open the NullGroup; it may fail if there is no matching
        # trace file. Skip if this is true.
        null = reconx.NullGroup(f, dir_swmf, imf)

        # Make sure we just have 1 null pair in our group.
        # if len(null) > 1:
        #     print(f'WARNING: MULTIPLE NULLS AT {null.time}')
        #     continue

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
        data['pram'].append(null.pram)
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


def potential_plots(runname):
    '''
    Given a results pickle, make a 3-panel plot of the varying SW driver, 
    the four potentials, and their correlations.

    Input: runname: string, one of the keys in the runs dictionary.
    '''

    # Load a pickle.
    with open(f'reconexus_results_{runname}.pkl', 'rb') as f:
        data = pickle.load(f)

    imf = ImfInput(runs[runname]['imffile'])
    hours_imf = np.array([(t - imf['time'][0]).total_seconds()/3600. 
                            for t in imf['time']])
    interp = np.interp(data['runtime'], hours_imf, imf['ux'])

    slope, intercept, r_value, p_value, std_err  = linregress(data['cpcp'], data['krpot'])

    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(16, 12))
    fig.suptitle(f'Results from {runname.capitalize()} Investigation')

    if runname == 'bz' or runname == 'bz_vsw750':
        a1.plot(data['runtime'], data['bz'], 'o-')
        a1.set_ylabel('IMF $B_z$')
    elif runname == 'vsw':
        a1.plot(data['runtime'], data['u'], 'o-')
        a1.set_ylabel('$V_{sw}$')
    elif runname == 'n':
        a1.plot(data['runtime'], data['n'], 'o-')
        a1.set_ylabel('$n_{sw}$')
    elif runname == 'cond':
        a1.plot(data['runtime'], data['ped'], 'o-')
        a1.set_ylabel('Uniform $\\Sigma_{Ped}$')
    else:
        raise ValueError(f"Unknown runname: {runname!r}. ")
    
    a1.set_xlabel(r'Hours From Simulation Start (S.T.)')
    a1.axvline(x=6, linestyle='--',color='black', label='Southward Turning')
    a1.legend(loc='best')

    a2.plot(data['runtime'], data['geopot'], color='Blue',
            marker='o', label='$\\Phi_{RXN}$')
    a2.plot(data['runtime'], data['cpcp'], color='Red',
            marker='o', label='CPCP')
    a2.plot(data['runtime'], data['krpot'], color='Gray',
            marker='o', label='$\\Phi_{K&R}$')
    a2.plot(data['runtime'], data['ifpot'], color='Green',
            marker='o', label='$\\Phi_{IFP}$')
    a2.axvline(x=6, linestyle='--',color='black')
    a2.legend(loc='best')
    a2.set_xlabel(r'Hours From Simulation Start (S.T.)')
    a2.set_ylabel(r'$\Phi$ ($kV$)')

    # CPCP vs DeltaV correlation
    a3.scatter(data['cpcp'], data['krpot'], color='Purple', label='Data')
    a3.plot(data['cpcp'],
            intercept + slope*np.array(data['cpcp']),
            color='Black', linestyle='--',
            label='Fit: $\\Phi_{K&R}$'
            +f'={slope:.2f}*CPCP+{intercept:.2f}\nR={r_value:.2f}')
    a3.set_xlabel('CPCP (SWMF) [kV]')
    a3.set_ylabel('$\\Phi_{K&R}$[kV]')
    a3.set_title('Correlation: CPCP vs $\\Phi_{K&R}$')
    a3.grid(True)
    a3.legend(loc='best', fontsize=8)
    
    fig.tight_layout()


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
    a2.plot(data['runtime'], data['alf'], 'o-', c='C1')
    a2.set_ylabel(r'$V_{Alf, SW}$', c='C1')

    a3.plot(data['runtime'], data['ped'], label=r'$\Sigma_{Ped,Iono}$')
    a3.plot(data['runtime'], data['sigmaA'], label=r'$\Sigma_{SW}$')
    a3.legend(loc='best')
    a3.set_ylabel('mhos')

    a4.plot(data['runtime'], data['gel'], 'o-')
    a4.set_ylabel('GEL ($R_E$)')

    a5.plot(data['runtime'], data['geopot'], color='Blue',
            marker='o', label=r'Reconexus $\Phi$')
    a5.plot(data['runtime'], data['cpcp'], color='Red',
            marker='o', label='CPCP')
    a5.plot(data['runtime'], data['krpot'], color='Gray',
            marker='o', label='K&R Pot')
    a5.plot(data['runtime'], data['ifpot'], color='Green',
            marker='o', label='Footpoint Pot')
    a5.legend(loc='best')
    a5.set_ylabel('$\\Phi$ ($kV$)')

    fig.tight_layout()


def viscous_potential_plot(runname):
    '''
    Given a results pickle, calculate the difference between CPCP and Iono Footpoint
    potentials to quantify viscous potential contribution. Plot this difference 
    against the varying SW drivers (namely vsw and n).
    '''

    # Open the files we want:
    # Load a pickle with reconnexus stuff:
    with open(f'reconexus_results_{runname}.pkl', 'rb') as f:
        data = pickle.load(f)
    
    viscpot = data['cpcp'] - np.array(data['ifpot'])

    if runname == 'vsw':
        xvar = data['u']
        xlabel = '$V_{sw}$ (km/s)'
        mask = (data['runtime'] >= 6) & (data['runtime'] <= 22)
        pram = data['pram'][mask]
        viscpot = viscpot[mask]
    elif runname == 'n':
        xvar = data['n']
        pram = data['pram']
        viscpot = viscpot
    else:
        raise ValueError(f"Unknown runname: {runname!r}. Must be 'vsw' or 'n'.")

    # Calculate linear regression for viscous potential vs RAM pressure
    slope, intercept, r_value, p_value, std_err = linregress(pram, viscpot)
    xfit = np.linspace(np.min(pram), np.max(pram), 100)
    yfit = slope*xfit + intercept

    fig, a1 = plt.subplots(figsize=(8,6))
    a1.scatter(pram, viscpot, marker='o', label='Simulation Data')

    a1.plot(xfit, yfit,  '--', color='Black', label=(f'Fit: $\Delta\Phi$ = {slope:.3f} {runname} + {intercept:.2f}\n'f'R = {r_value:.2f}'))
    a1.set_xlabel('$P_{dyn}$ (nPa)', fontsize=16)
    a1.set_ylabel('$\\Phi_{CPCP}-\\Phi_{IFP}$ (kV)', fontsize=16)
    a1.set_title('Excess Potential vs Dynamic Pressure', fontsize=18)

    a1.grid(True)
    a1.legend()

    plt.tight_layout()


def compare_gel_2times(runname, t1=dt.datetime(1998, 5, 4, 0, 0, 0),
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
    dirrxn = runs[runname]['rxn_dir']
    dirswmf = runs[runname]['swmf_dir']

    # Open our RxN lines, separators, and MHD files:
    strtime = f"{t1:%Y%m%d-%H%M%S}"
    lp1 = reconx.read_nulls(dirrxn + f"null_line_pls_n01_001_e{strtime}.dat")
    ln1 = reconx.read_nulls(dirrxn + f"null_line_neg_n01_001_e{strtime}.dat")
    sep1 = reconx.read_separator(dirrxn + f"Separator_e{strtime}.dat")
    mhd1 = bats.Bats2d(dirswmf + f"GM/z=0_mhd_2_e{strtime}-000.out")
    mhd1.calc_j()

    strtime = f"{t2:%Y%m%d-%H%M%S}"
    lp2 = reconx.read_nulls(dirrxn + f"null_line_pls_n01_001_e{strtime}.dat")
    ln2 = reconx.read_nulls(dirrxn + f"null_line_neg_n01_001_e{strtime}.dat")
    sep2 = reconx.read_separator(dirrxn + f"Separator_e{strtime}.dat")
    mhd2 = bats.Bats2d(dirswmf + f"GM/z=0_mhd_2_e{strtime}-000.out")

    kwargs = {'ylim': [-45, 15], 'xlim': [-25, 25], 'cmap': 'viridis',
              'nlev': 51, 'dolog': True, 'extend': 'max', 'add_cbar': True}


    fig = plt.figure(figsize=[16, 12])
    a1, a2 = fig.add_subplot(2, 2, 1), fig.add_subplot(2, 2, 2, projection='3d')
    a3 = fig.add_subplot(2, 1, 2)
    mhd1.add_contour('y', 'x', 'j', target=a1, loc=121, **kwargs)
    a1.plot(lp1['Y'], lp1['X'], '.r', ln1['Y'], ln1['X'], 
            'ro', label=f'Time 1 = {t1:%d-%H%M}')
    a1.plot(lp2['Y'], lp2['X'], '.b', ln2['Y'], ln2['X'], 
            'bo', label=f'Time 2 = {t2:%d-%H%M}')
    a1.plot(sep1['Y'], sep1['X'], 'mx', label='Separator 1')
    a1.plot(sep2['Y'], sep2['X'], 'gx', label='Separator 2')
    a1.legend(loc='best')

    a2.plot(lp1['X'], lp1['Y'], lp1['Z'], '.r')
    a2.plot(ln1['X'], ln1['Y'], ln1['Z'], '.r', label='Time 1')
    a2.plot(lp2['X'], lp2['Y'], lp2['Z'], '.b')
    a2.plot(ln2['X'], ln2['Y'], ln2['Z'], '.b', label='Time 2')
    a2.plot(sep1['X'], sep1['Y'], sep1['Z'], 'mx', label='Separator 1')
    a2.plot(sep2['X'], sep2['Y'], sep2['Z'], 'gx', label='Separator 2')

    a3.plot(data['runtime'], data['geopot'], color='Blue', marker='o', 
            linestyle='-', label=r'$\Phi_{RECONEXUS}$')
    a3.plot(data['runtime'], data['cpcp'], color='Red', marker='o', 
            linestyle='-', label='CPCP')
    a3.plot(data['runtime'], data['krpot'], color='Gray', marker='o', 
            linestyle='-', label=r'$\Phi_{K&R}$')
    a3.plot(data['runtime'], data['ifpot'], color='Green', marker='o', 
            linestyle='-', label=r'$\Phi_{IFP}$')
    a3.set_xlabel(r'Hours from Simulation Start (S.T.)')
    a3.set_ylabel(r'$\Phi$ ($kV$)')
    a3.legend(loc='best')


def compare_2runs_2times(runname1, runname2, t1=dt.datetime(1998, 5, 4, 8, 0, 0),
                       t2=dt.datetime(1998, 5, 5, 11, 0, 0)):
    '''
    Examine the position of separators for two runs at two times.

    t1 and t2 should be datetimes. :P
    '''
    # Open the files we want:
    # Load a pickle with reconnexus stuff:
    with open(f'reconexus_results_{runname1}.pkl', 'rb') as f:
        data1 = pickle.load(f)

    with open(f'reconexus_results_{runname2}.pkl', 'rb') as f:
        data2 = pickle.load(f)

    # Shortcut vars:
    dirrxn1, dirswmf1 = runs[runname1]['rxn_dir'], runs[runname1]['swmf_dir']
    dirrxn2, dirswmf2 = runs[runname2]['rxn_dir'], runs[runname2]['swmf_dir']

    # Open our RxN lines, separators, and MHD files:
    strtime = f"{t1:%Y%m%d-%H%M%S}"
    lpt1r1 = reconx.read_nulls(dirrxn1 + f"null_line_pls_n01_001_e{strtime}.dat")
    lnt1r1 = reconx.read_nulls(dirrxn1 + f"null_line_neg_n01_001_e{strtime}.dat")
    sept1r1 = reconx.read_separator(dirrxn1 + f"Separator_e{strtime}.dat")
    mhdt1r1 = bats.Bats2d(dirswmf1 + f"GM/z=0_mhd_2_e{strtime}-000.out")
    mhdt1r1.calc_j()

    lpt1r2 = reconx.read_nulls(dirrxn2 + f"null_line_pls_n01_001_e{strtime}.dat")
    lnt1r2 = reconx.read_nulls(dirrxn2 + f"null_line_neg_n01_001_e{strtime}.dat")
    sept1r2 = reconx.read_separator(dirrxn2 + f"Separator_e{strtime}.dat")
    mhdt1r2 = bats.Bats2d(dirswmf2 + f"GM/z=0_mhd_2_e{strtime}-000.out")
    mhdt1r2.calc_j()

    strtime = f"{t2:%Y%m%d-%H%M%S}"
    lpt2r1 = reconx.read_nulls(dirrxn1 + f"null_line_pls_n01_001_e{strtime}.dat")
    lnt2r1 = reconx.read_nulls(dirrxn1 + f"null_line_neg_n01_001_e{strtime}.dat")
    sept2r1 = reconx.read_separator(dirrxn1 + f"Separator_e{strtime}.dat")
    mhdt2r1 = bats.Bats2d(dirswmf1 + f"GM/z=0_mhd_2_e{strtime}-000.out")
    mhdt2r1.calc_j()

    lpt2r2 = reconx.read_nulls(dirrxn2 + f"null_line_pls_n01_001_e{strtime}.dat")
    lnt2r2 = reconx.read_nulls(dirrxn2 + f"null_line_neg_n01_001_e{strtime}.dat")
    sept2r2 = reconx.read_separator(dirrxn2 + f"Separator_e{strtime}.dat")
    mhdt2r2 = bats.Bats2d(dirswmf2 + f"GM/z=0_mhd_2_e{strtime}-000.out")
    mhdt2r2.calc_j()

    kwargs = {'ylim': [-45, 15], 'xlim': [-25, 25], 'cmap': 'viridis', 'nlev': 51,
        'dolog': True, 'extend': 'max', 'add_cbar': True}

    fig = plt.figure(figsize=[16, 16])

    a1, a2 = fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2, projection='3d')
    mhdt1r1.add_contour('y', 'x', 'j', target=a1, loc=121, **kwargs)
    mhdt1r2.add_contour('y', 'x', 'j', target=a1, loc=121, **kwargs, alpha=0.5)
    a1.plot(lpt1r1['Y'], lpt1r1['X'], '^r', label=f'{runname1}, {t1:%d-%H%M}')
    a1.plot(lnt1r1['Y'], lnt1r1['X'], '^r') 
    a1.plot(lpt2r1['Y'], lpt2r1['X'], '^b', label=f'{runname1}, {t2:%d-%H%M}')
    a1.plot(lnt2r1['Y'], lnt2r1['X'], '^b')
    a1.plot(sept1r1['Y'], sept1r1['X'], '^m', label=f'{runname1}, {t1:%d-%H%M}')
    a1.plot(sept2r1['Y'], sept2r1['X'], '^y', label=f'{runname1}, {t2:%d-%H%M}')
    a1.plot(lpt1r2['Y'], lpt1r2['X'], '.r', label=f'{runname2}, {t1:%d-%H%M}')
    a1.plot(lnt1r2['Y'], lnt1r2['X'], '.r')
    a1.plot(lpt2r2['Y'], lpt2r2['X'], '.b', label=f'{runname2}, {t2:%d-%H%M}')
    a1.plot(lnt2r2['Y'], lnt2r2['X'], '.b')
    a1.plot(sept1r2['Y'], sept1r2['X'], '.m',  label=f'{runname2}, {t1:%d-%H%M}')
    a1.plot(sept2r2['Y'], sept2r2['X'], '.y',  label=f'{runname2}, {t2:%d-%H%M}')
    a1.legend(loc='best')

    a2.plot(lpt1r1['X'], lpt1r1['Y'], lpt1r1['Z'], '.r')
    a2.plot(lnt1r1['X'], lnt1r1['Y'], lnt1r1['Z'], '.r', 
            label=f'{runname1} {t1:%d-%H%M}')
    a2.plot(lpt2r1['X'], lpt2r1['Y'], lpt2r1['Z'], '.b')
    a2.plot(lnt2r1['X'], lnt2r1['Y'], lnt2r1['Z'], '.b', 
            label=f'{runname1} {t2:%d-%H%M}')
    a2.plot(sept1r1['X'], sept1r1['Y'], sept1r1['Z'], '.m', 
            label=f'{runname1},{t1:%d-%H%M}')
    a2.plot(sept2r1['X'], sept2r1['Y'], sept2r1['Z'], '.y', 
            label=f'{runname1}, {t2:%d-%H%M}')
    a2.plot(lpt1r2['X'], lpt1r2['Y'], lpt1r2['Z'], '^r')
    a2.plot(lnt1r2['X'], lnt1r2['Y'], lnt1r2['Z'], '^r', 
            label=f'{runname2} {t1:%d-%H%M}')
    a2.plot(lpt2r2['X'], lpt2r2['Y'], lpt2r2['Z'], '^b')
    a2.plot(lnt2r2['X'], lnt2r2['Y'], lnt2r2['Z'], '^b', 
            label=f'{runname2} {t2:%d-%H%M}')
    a2.plot(sept1r2['X'], sept1r2['Y'], sept1r2['Z'], '^m', 
            label=f'{runname2}, {t1:%d-%H%M}')
    a2.plot(sept2r2['X'], sept2r2['Y'], sept2r2['Z'], '^y', 
            label=f'{runname2}, {t2:%d-%H%M}')

def plot_nulls_and_separator(runname, t=dt.datetime(1998, 5, 4, 0, 0, 0)):
    '''
    Given a runname, retrieve the null and separator file to plot the nulls and 
    endpoints in 3D
    '''

    dirrxn, dirswmf = runs[runname]['rxn_dir'], runs[runname]['swmf_dir']

    strtime = f"{t:%Y%m%d-%H%M%S}"

    mhdz = bats.Bats2d(dirswmf + f"GM/z=0_mhd_2_e{strtime}-000.out")
    mhdz.calc_j()

    # Open our RxN lines, separators, and MHD files:
    strtime = f"{t:%Y%m%d-%H%M%S}"
    pnulls = reconx.read_nulls(dirrxn + f"PlusNulls_e{strtime}.dat")
    nnulls = reconx.read_nulls(dirrxn + f"NegNulls_e{strtime}.dat")

    #sep = reconx.read_separator(dirrxn + f"Separator_e{strtime}.dat")
    ptrace = reconx.read_nulls(dirrxn + f"null_line_pls_n01_001_e{strtime}.dat", reorder=True)
    ntrace = reconx.read_nulls(dirrxn + f"null_line_neg_n01_001_e{strtime}.dat", reorder=True)

    fig = plt.figure(figsize=[16, 16])
    a1, a2 = fig.add_subplot(1, 2, 1, projection='3d'), fig.add_subplot(1,2,2)
    mhdz.add_contour('y', 'x', 'j', dolog=True, target=a2, loc=112, ylim=[-60,20])
    a2.plot(pnulls['Y'], pnulls['X'], 'xr')
    a2.plot(nnulls['Y'], nnulls['X'], 'xb')
    a1.plot(0, 0, 0, 'ok')
    a1.plot(ptrace['X'], ptrace['Y'], ptrace['Z'], '.r')
    a1.plot(ntrace['X'], ntrace['Y'], ntrace['Z'], '.b')
    a1.plot(pnulls['X'], pnulls['Y'], pnulls['Z'], 'xr')
    a1.plot(nnulls['X'], nnulls['Y'], nnulls['Z'], 'xb')
    a1.set_aspect('equal')
    #a1.plot(sep['X'], sep['Y'], sep['Z'], '.b')
    a1.set_xlabel('X')
    a1.set_ylabel('Y')
    a1.set_zlabel('Z')
    a1.set_title(f"At time= {t}")
    fig.tight_layout()

        #plt.savefig()