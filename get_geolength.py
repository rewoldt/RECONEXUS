#!/usr/bin/env python                                                                                                       

'''                                                                                                                          
Some simple tools to fit ray traces and calculate geoeffective length                                                        
'''

import numpy as np
import matplotlib.pyplot as plt
import reconx
import get_footprints
import get_iono_drop
import os
import re
from spacepy.pybats import bats, rim
from datetime import datetime
from scipy.stats import pearsonr
import glob

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

def get_potential(null_file_neg, null_file_pls, n):
    '''
    Calculate geoeffective length and potential given a set of field lines.

    Parameters
    ----------
    null_file_neg : Field lines starting from negative null points.
    null_file_pls : Field lines starting from positive null points.
    n : The index to slice from.
    u : The velocity field vector.
    B : The magnetic field vector.

    Returns
    -------
    dist : List of distances between nulls.
    geo_potential : Geoeffective potential calculated as the dot product of u cross B and the null-to-null vector.
    '''
    
    Xp_sorted, Yp_sorted, Zp_sorted = get_footprints.sort_trace(null_file_pls)
    Xn_sorted, Yn_sorted, Zn_sorted = get_footprints.sort_trace(null_file_neg)
    
    # Slicing the field lines from index n onwards
    Xn = Xn_sorted[n:]
    Yn = Yn_sorted[n:]
    Zn = Zn_sorted[n:]
    Xp = Xp_sorted[n:]
    Yp = Yp_sorted[n:]
    Zp = Zp_sorted[n:]
    
    RE2km = 6378.14  # Radius of Earth in km
    
    dist = []
    length = []
    
    # Calculate u x B (cross product of velocity and magnetic field)
    uxb = np.cross(u, B)
    
    for i in range(len(Xn)):
        for j in range(len(Xp)):
            # Calculate the distance vector between corresponding nulls
            dist_vector = [(Xp[i] - Xn[j]) * RE2km, 
                           (Yp[i] - Yn[j]) * RE2km, 
                           (Zp[i] - Zn[j]) * RE2km]
            
            # Calculate the magnitude of the distance vector
            distance = np.linalg.norm(dist_vector)
            dist.append(distance)
            length.append(dist_vector)
    
    # Find the index of the minimum distance
    min_distance_index = np.argmin(dist)
    min_distance_vector = length[min_distance_index]
    
    # Calculate the geoeffective potential as the dot product of u x B and the minimum distance vector
    geo_potential = np.dot(uxb, min_distance_vector)
    
    return dist, abs(geo_potential)

def plot_tracers_3D(trace_file_neg, trace_file_pls, null_file_neg, null_file_pls, null_number):
    # Read the tracer points from the files
    
    # Extract the null points based on null number
    Xn_tracer = trace_file_neg['X']
    Yn_tracer = trace_file_neg['Y']
    Zn_tracer = trace_file_neg['Z']
    Xp_tracer = trace_file_pls['X']
    Yp_tracer = trace_file_pls['Y']
    Zp_tracer = trace_file_pls['Z']
    

    # Read the null points from the files
    null_neg = reconx.read_nulls(null_file_neg)
    null_pls = reconx.read_nulls(null_file_pls)

    # Extract the null points based on null number
    Xn_null = null_neg['X'][null_number]
    Yn_null = null_neg['Y'][null_number]
    Zn_null = null_neg['Z'][null_number]
    Xp_null = null_pls['X'][null_number]
    Yp_null = null_pls['Y'][null_number]
    Zp_null = null_pls['Z'][null_number]

    # Plotting in 3D
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # Origin point for reference
    ax.scatter3D(0, 0, 0, s=50, color='Black')
    # Plot traces
    ax.scatter3D(Xn_tracer, Yn_tracer, Zn_tracer, color='Blue', label='Negative trace')
    ax.scatter3D(Xp_tracer, Yp_tracer, Zp_tracer, color='Red', label='Positive trace')
    # Plot null points
    ax.scatter3D(Xn_null, Yn_null, Zn_null, color='Black', marker='x', s=50, label='Negative null')
    ax.scatter3D(Xp_null, Yp_null, Zp_null, color='Black', marker='x', s=50, label='Positive null')

    # Setting labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.legend()
    #plt.text()
    #plt.show()
    plt.savefig('./plots/AMR_40000_plots/3D_AMR_'+str(null_number)+'.png', dpi=300)

def plot_tracers_2D(trace_file_neg, trace_file_pls, null_file_neg, null_file_pls, null_number):
    # Read the tracer points from the files
    
    # Extract the null points based on null number
    Xn_tracer = trace_file_neg['X']
    Yn_tracer = trace_file_neg['Y']
    Zn_tracer = trace_file_neg['Z']
    Xp_tracer = trace_file_pls['X']
    Yp_tracer = trace_file_pls['Y']
    Zp_tracer = trace_file_pls['Z']
    

    # Read the null points from the files
    null_neg = reconx.read_nulls(null_file_neg)
    null_pls = reconx.read_nulls(null_file_pls)

    # Extract the null points based on null number
    Xn_null = null_neg['X'][null_number]
    Yn_null = null_neg['Y'][null_number]
    Zn_null = null_neg['Z'][null_number]
    Xp_null = null_pls['X'][null_number]
    Yp_null = null_pls['Y'][null_number]
    Zp_null = null_pls['Z'][null_number]
    print (Xn_null, Yn_null, Zn_null, Xp_null, Yp_null, Zp_null)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))  # Create two subplots side by side
    
    # First subplot: y=0
    ax1.plot(0, 0, marker='o', markersize=5, color='Black')  # Change 's' to 'marker' and 'markersize' for plot
    # Plot traces
    ax1.scatter(Xn_tracer, Zn_tracer, color='Blue', label='Negative trace')
    ax1.scatter(Xp_tracer, Zp_tracer, color='Red', label='Positive trace')
    # Plot null points
    ax1.scatter(Xn_null, Zn_null, color='Black', marker='x', s=50, label='Negative null')
    ax1.scatter(Xp_null, Zp_null, color='Black', marker='x', s=50, label='Positive null')
    ax1.set_title('y=0')
    ax1.set(xlabel='X (RE)', ylabel='Z (RE)')
    ax1.set_aspect('equal')
    
    # Second subplot: z=0
    ax2.plot(0, 0, marker='o', markersize=5, color='Black')  # Change 's' to 'marker' and 'markersize' for plot
    # Plot traces
    ax2.scatter(Xn_tracer, Yn_tracer, color='Blue', label='Negative trace')
    ax2.scatter(Xp_tracer, Yp_tracer, color='Red', label='Positive trace')
    # Plot null points
    ax2.scatter(Xn_null, Yn_null, color='Black', marker='x', s=50, label='Negative null')
    ax2.scatter(Xp_null, Yp_null, color='Black', marker='x', s=50, label='Positive null')
    ax2.set_title('z=0')
    ax2.set(xlabel='X (RE)', ylabel='Y (RE)')
    ax2.set_aspect('equal')
    
    # Display legends and show plots
    ax1.legend()
    plt.tight_layout()  # Adjust the layout to prevent overlap
    #plt.show()
    plt.savefig('./plots/AMR_40000_plots/2D_AMR_'+str(null_number)+'.png', dpi=300)

def plot_paneled_geo_potential(trace_data_neg, trace_data_pls, ionosphere_data, null_file_neg, null_file_pls, null_number, timestamp):
    # Read the null points from the files
    null_neg = reconx.read_nulls(null_file_neg)
    null_pls = reconx.read_nulls(null_file_pls)

    # Extract the null points based on null number
    Xn_null = null_neg['X'][null_number]
    Yn_null = null_neg['Y'][null_number]
    Zn_null = null_neg['Z'][null_number]
    Xp_null = null_pls['X'][null_number]
    Yp_null = null_pls['Y'][null_number]
    Zp_null = null_pls['Z'][null_number]
    
    # Create a paneled plot with 3 subplots: 2D cuts and a polar plot
    fig = plt.figure(figsize=(15, 5))
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133, projection='polar')

    # First subplot: y=0 2D cut
    ax1.scatter(0, 0, color='Black')
    ax1.scatter(trace_data_neg['X'], trace_data_neg['Z'], color='Blue', label='Negative trace')
    ax1.scatter(trace_data_pls['X'], trace_data_pls['Z'], color='Red', label='Positive trace')
    ax1.scatter(Xn_null, Zn_null, color='Black', marker='x', s=50, label='Negative null')
    ax1.scatter(Xp_null, Zp_null, color='Black', marker='x', s=50, label='Positive null')
    
    ax1.set_title('y=0')
    ax1.set(xlabel='X (RE)', ylabel='Z (RE)')
    ax1.set_aspect('equal')
    ax1.legend()

    # Add text at the bottom-left corner of the first plot
    #ax1.text(0, 0, 'Text for ax1', transform=ax1.transAxes, ha='left', va='bottom')

    # Second subplot: z=0 2D cut
    
    ax2.scatter(0, 0, color='Black')
    ax2.scatter(trace_data_neg['X'], trace_data_neg['Y'], color='Blue', label='Negative trace')
    ax2.scatter(trace_data_pls['X'], trace_data_pls['Y'], color='Red', label='Positive trace')
    ax2.scatter(Xn_null, Yn_null, color='Black', marker='x', s=50, label='Negative null')
    ax2.scatter(Xp_null, Yp_null, color='Black', marker='x', s=50, label='Positive null')
    ax2.set_title('z=0')
    ax2.set(xlabel='X (RE)', ylabel='Y (RE)')
    ax2.set_aspect('equal')
    #ax2.legend()

    # Add calculated potential at the bottom-left corner of the second plot
    dist, geo_potential = get_potential(trace_data_neg, trace_data_pls, -50)
    ax2.text(0.5, 0.05, f'GEL Potential={geo_potential:.2f}$kV$', transform=ax2.transAxes)

    # Third subplot: Ionosphere contour plot
    ionosphere_data.add_cont('n_phi', target=ax3, add_cbar=True)  # Use ax3 for the contour plot
    ax3.set_title('Ionosphere Potential Contour')

    # Add calculated potentials (CPCP and potential drop) at the bottom-left corner of the third plot
    Xp_sorted, Yp_sorted, Zp_sorted = get_footprints.sort_trace(trace_data_pls)
    Xn_sorted, Yn_sorted, Zn_sorted = get_footprints.sort_trace(trace_data_neg)
    coord_pls = get_footprints.get_footprints(Xp_sorted, Yp_sorted, Zp_sorted)
    coord_neg = get_footprints.get_footprints(Xn_sorted, Yn_sorted, Zn_sorted)
    pot_results = get_iono_drop.get_iono_drop(ionosphere_data, coord_neg, coord_pls)

    ax3.text(-0.1, -0.05, f'CPCP={pot_results["cpcp"]:.2f}$kV$\nNull-Null Drop={pot_results["potdrop"]:.2f}$kV$', transform=ax3.transAxes)

    plt.tight_layout()
    plt.savefig('./plots_psphere/paneled_geo_potential_'+str(timestamp)+'_'+str(null_number)+'.png', dpi=300)
    #plt.show()

    
def plot_current(mhd_file, null_file_neg, null_file_pls, timestamp):

    '''
    
    Plots nulls on top of mhd current plot
    
    '''
    # Read the null points from the files
    null_neg = reconx.read_nulls(null_file_neg)
    null_pls = reconx.read_nulls(null_file_pls)

    # Extract the null points based on null number
    Xn_null = null_neg['X']
    Yn_null = null_neg['Y']
    Zn_null = null_neg['Z']
    Xp_null = null_pls['X']
    Yp_null = null_pls['Y']
    Zp_null = null_pls['Z']
    
    mhd = bats.Bats2d(mhd_file)
    t = mhd.attrs['strtime']
    print (t)
    mhd.calc_j()
    fig = plt.figure(figsize=[10,10])
    fig = mhd.add_contour('x', 'y', 'j', target=fig, dolog=True, xlim=[-60, 20], ylim=[-40, 40])
    plt.scatter(Xn_null, Yn_null, color='Blue', marker='x', label='Negative null')
    plt.scatter(Xp_null, Yp_null, color='Red', marker='x', label='Positive null')
    plt.title(f'time={t}')
    
    # Connect corresponding positive and negative nulls with lines
    for x_neg, y_neg, x_pos, y_pos in zip(Xn_null, Yn_null, Xp_null, Yp_null):
        plt.plot([x_neg, x_pos], [y_neg, y_pos], color='green', linestyle='--', linewidth=1)

    plt.tight_layout()
    plt.savefig('./plots_psphere/test_nulls_'+str(timestamp)+'.png', dpi=300)

def plot_time_dynamics(results):
    # Sort the timestamps
    sorted_timestamps = sorted(results.keys(), key=lambda x: datetime.strptime(x, '%H%M%S'))
    
    plt.figure()
    plt.style.use('seaborn-v0_8')
    
    # Iterate over each null number
    for null_number in results[sorted_timestamps[0]].keys():
        cpcp_values = []
        potdrop_values = []
        geo_potential_values = []
        times = []
        
        for timestamp in sorted_timestamps:
            if null_number in results[timestamp]:
                try:
                    times.append(datetime.strptime(timestamp, '%H%M%S'))
                    cpcp_values.append(np.mean([results[timestamp][null_number][null_tag]['cpcp'] for null_tag in results[timestamp][null_number]]))
                    potdrop_values.append(np.mean([results[timestamp][null_number][null_tag]['potdrop'] for null_tag in results[timestamp][null_number]]))
                    geo_potential_values.append(np.mean([results[timestamp][null_number][null_tag]['geo_potential'] for null_tag in results[timestamp][null_number]]))
                except IndexError:
                    # Handle missing data
                    print(f"Missing data for {null_number} at {timestamp}")
                    continue
        
        # Plotting
        plt.plot(times, cpcp_values, label="CPCP", marker='*', linestyle='-', color='b')
        #plt.scatter(times, potdrop_values, label="Potential Drop Null", marker='*', linestyle='--', color='r')
        plt.plot(times, geo_potential_values, label="GEL Potential", marker='^', linestyle='-', color='r')

    plt.xlabel('Time')
    plt.ylabel('Potential (kV)')
    plt.title('Potential vs Time')
    plt.grid(True)
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Moving the legend to the top right and making it multiple columns
    plt.legend()
    
    # Save the plot if needed
    plt.savefig('./plots_psphere/time_dynamics_psphere.png', dpi=300)
    
    # Show the plot
    #plt.show()

'''
def plot_results(results):
    # Sort the timestamps
    sorted_timestamps = sorted(results.keys(), key=lambda x: datetime.strptime(x, '%H%M%S'))
    
    # Collect all null numbers
    unique_null_numbers = {null_number for timestamp_data in results.values() for null_number in timestamp_data.keys()}
    
    for null_number in unique_null_numbers:
        times = []
        mean_cpcp_values = []
        mean_potdrop_values = []
        mean_geo_potential_values = []

        for timestamp in sorted_timestamps:
            if null_number in results[timestamp]:
                try:
                    times.append(timestamp)
                    
                    # Aggregate data for all null tags under this null number and timestamp
                    all_cpcp = []
                    all_potdrop = []
                    all_geo_potential = []
                    
                    for null_tag in results[timestamp][null_number]:
                        all_cpcp.extend(results[timestamp][null_number][null_tag]['cpcp'])
                        all_potdrop.extend(results[timestamp][null_number][null_tag]['potdrop'])
                        all_geo_potential.extend(results[timestamp][null_number][null_tag]['geo_potential'])
                    
                    # Calculate the mean values
                    mean_cpcp_values.append(np.mean(all_cpcp))
                    mean_potdrop_values.append(np.mean(all_potdrop))
                    mean_geo_potential_values.append(np.mean(all_geo_potential))
                except IndexError:
                    # Handle missing data
                    print(f"Missing data for {null_number} at {timestamp}")
                    continue
        
        # Plotting
        #plt.semilogy(times, mean_cpcp_values, marker='^', linestyle=':', color='b')
        #plt.semilogy(times, mean_potdrop_values, marker='s', linestyle='--', color='g')
        #plt.semilogy(times, mean_geo_potential_values, marker='o', color='r')
    plt.figure(figsize=(12, 8))
    #plt.scatter(mean_geo_potential_values, mean_cpcp_values, label=f"Null Number {null_number}", marker='^', color='b')

        # Calculate correlation coefficient
    #correlation, _ = pearsonr(mean_geo_potential_values, mean_cpcp_values)
    
    # Add text for correlation
    #plt.text(0.05, 0.95, f'Correlation: {correlation:.2f}', ha='left', va='center', transform=plt.gca().transAxes)
    plt.plot(times, mean_cpcp_values, label=f'CPCP', marker='^', linestyle=':', color='b')
    plt.plot(times, mean_potdrop_values, label=f'Potential drop', marker='s', linestyle='--', color='g')
    plt.plot(times, mean_geo_potential_values, label=f'GEL Potential', marker='o', color='r')
    # Customize the plot
    plt.xlabel('GEL Potential')
    plt.ylabel('CPCP (kV)')
    plt.title('Comparison of  CPCP and GEL Potential for All Null Pairs')
    plt.legend(loc='best', ncol=3)#, prop={'size': 5.5})
    plt.grid(True)
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save the plot if needed
    plt.savefig('./time_dynamics_amr.png', dpi=300)
'''

def plot_mean_cpcp_vs_geo_potential(results):
    mean_cpcp_values = []
    mean_geo_potential_values = []
    mean_potdrop_values = []
    times = []
    null_pairs = []

    for timestamp in results:
        times.append(datetime.strptime(timestamp, '%H%M%S'))
        for null_number in results[timestamp]:
            cpcp_values = []
            geo_potential_values = []
            potdrop_values = []
            
            
            for null_tag in results[timestamp][null_number]:
                cpcp_values.extend(results[timestamp][null_number][null_tag]['cpcp'])
                geo_potential_values.extend(results[timestamp][null_number][null_tag]['geo_potential'])
                potdrop_values.extend(results[timestamp][null_number][null_tag]['potdrop'])

            if cpcp_values and geo_potential_values and potdrop_values:
                mean_cpcp = np.mean(cpcp_values)
                mean_geo_potential = np.mean(geo_potential_values)
                mean_potdrop = np.mean(potdrop_values)
                
                mean_cpcp_values.append(mean_cpcp)
                mean_geo_potential_values.append(mean_geo_potential)
                mean_potdrop_values.append(mean_potdrop)
                null_pairs.append(null_number)

    # Plotting
    plt.figure(figsize=(12, 8))
    plt.style.use('seaborn-v0_8')
    
    correlation_cpcp, _ = pearsonr(mean_geo_potential_values, mean_cpcp_values)
    # Plotting Mean CPCP vs Mean Geo Potential
    plt.scatter(mean_geo_potential_values, mean_cpcp_values, marker='.', color='b', label=f'Correlation: {correlation_cpcp:.2f}')
    #plt.scatter(mean_geo_potential_values, mean_cpcp_values, s=70, label='CPCP vs GEL Potential')
    # Calculate correlation coefficient for Mean CPCP vs Mean Geo Potential
    # Add text for correlation for Mean CPCP vs Mean Geo Potential
    #plt.text(0.05, 0.95, f'CPCP Correlation: {correlation_cpcp:.2f}', transform=plt.gca().transAxes)

    # Plotting Mean Potdrop vs Mean Geo Potential
    #plt.scatter(mean_geo_potential_values, mean_potdrop_values, marker='8', color='r', s=80, label='Iono Footprint Potential Drop vs GEL Potential')
    #plt.scatter(mean_geo_potential_values, mean_potdrop_values, s=70, color='indianred', label='Iono Footprint Potential Drop vs GEL Potential')
    #plt.scatter(times, mean_cpcp_values, label=f'CPCP', marker='^', linestyle='--', color='b')
    #plt.scatter(times, mean_potdrop_values, label=f'Potential drop', marker='s', linestyle='--', color='g')
    #plt.scatter(times, mean_geo_potential_values, label=f'GEL Potential', marker='o', linestyle='--', color='r')
    # Calculate correlation coefficient for Mean Potdrop vs Mean Geo Potential
    #correlation_potdrop, _ = pearsonr(mean_geo_potential_values, mean_potdrop_values)
    # Add text for correlation for Mean Potdrop vs Mean Geo Potential
    #plt.text(0.05, 0.90, f'Potdrop Correlation: {correlation_potdrop:.2f}', ha='left', va='center', transform=plt.gca().transAxes)

    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('GEL Potential (kV)', fontsize=16)
    plt.ylabel('CPCP (kV)', fontsize=16)
    plt.title('CPCP vs. GEL Potential', fontsize=18)
    plt.grid(True)
    plt.legend(fontsize=14)
    plt.tight_layout()

    # Save the plot
    plt.savefig('./plots_psphere/potential_comparison.png', dpi=300)

    # Show the plot
    #plt.show()

if __name__ == '__main__':
    # Your initial code goes here to populate the results dictionary
    # ux = -400 * 10**3  # m/s
    # uy = 0
    # uz = 0
    
    # Bx = 0
    # By = -2 * 10**-9  # T
    # Bz = -10 * 10**-9  # T
    
    ux = -450 * 10**3  # m/s
    uy = 0
    uz = 0
    
    Bx = 0
    By = 2 * 10**-9  # T
    Bz = -10 * 10**-9  # T
    
    B = [Bx, By, Bz]
    u = [ux, uy, uz]

    # Open a file:
    #iono = rim.Iono("/Users/rewoldt/HighResAMRdemo-ReconxWork/IE/it150321_000100_000.idl")
    #Directories for for the psphere runs
    #directory_ie = "../GmPlasSphere/plas_psgmie_smallby/"
    directory_ie = "/Volumes/coe-dwelling/swmf_results/GmPlasSphere/plas_psgmie_smallby/"
    directory = "/Volumes/coe-dwelling/swmf_results/reconnection_perfection/RECONX/run/"
    pattern = re.compile(r"null_line_neg_o\d+_\d{3}_t\d{8}\.dat")
    
    # Define the directory containing the trace files
    #directory = "/Users/rewoldt/RECONX/run_amr/"
    
    #Directory for locally store AMR runs
    #directory_ie = "/Users/rewoldt/HighResAMRdemo-ReconxWork/IE/"
    #directory = "/Users/rewoldt/RECONX/run/"
    #pattern = re.compile(r"null_line_neg_o\d+_\d{3}_e\d+-\d{6}\.dat")
    #For p'sphere runs
    # directory = "/Users/rewoldt/RECONX/run-psphere/"
    # pattern = re.compile(r"null_line_neg_o\d+_\d{3}_n\d+.dat") #null_line_neg_o03_001_n03030603.dat
    # iono = rim.Iono("./IE/it980504_160000_000.idl")
    #pattern = re.compile(r"null_line_neg_[ns]\d+_\d{3}_e\d+-\d{6}\.dat")
    #pattern = re.compile(r"null_line_neg_[ns]\d+_001_e\d+-\d{6}\.dat")
    results = {}
    
    # Loop over each trace file in the directory
    for filename in os.listdir(directory):
        #print (filename)
        if filename.endswith(".dat") and pattern.match(filename):
            # Extract timestamp from the filename
            #timestamp = filename.split("_")[-1].split("-")[-1].split(".")[0]
            #timestamp for psphere runs
            timestamp = filename.split("_")[-1][3:].split(".")[0]
            
            # Extract null number and polarity from the filename
            null_number = filename.split("_")[3][1:]
            polarity = filename.split("_")[3][:1]
            null_tag = filename.split("_")[4]
            #For AMR runs
            #nulls_neg = directory+f"NegNulls_e20150321-{timestamp}.dat"
            #nulls_pls = directory+f"PlusNulls_e20150321-{timestamp}.dat"
            #For psphere runs
            nulls_neg = directory+f'NegNulls_t00{timestamp}.dat'
            nulls_pls = directory+f'PlusNulls_t00{timestamp}.dat'
            #iono = rim.Iono(directory_ie+f'it150321_{timestamp}_000.idl')
            if timestamp == '000000': continue
            iono = rim.Iono(directory_ie+f'IE/it980504_{timestamp}_000.idl.gz')
            #nulls_neg = directory+f"NegNulls_n03030603.dat"
            #nulls_pls = directory+f"PlusNulls_n03030603.dat"
            null_number_plot = int(null_number) - 1
            
            try:
                # Read negative and positive trace files
                ray_neg = reconx.read_nulls(os.path.join(directory, filename))
                ray_pls = reconx.read_nulls(os.path.join(directory, filename.replace("neg", "pls")))
                
                # Calculate potential drop for the positive/negative pair
                dist, geo_potential = get_potential(ray_neg, ray_pls, -50)
                #dist, geo_potential = get_potential(ray_neg, ray_pls, 1)
    
                # Initialize the dictionary for this timestamp if it doesn't exist
                if timestamp not in results:
                    results[timestamp] = {}
                
                # Initialize the dictionary for this null_number if it doesn't exist
                if null_number not in results[timestamp]:
                    #results[timestamp][null_number] = {"geo_potential": [], "cpcp": [], "potdrop": []}
                    results[timestamp][null_number] = {}
                    #results[timestamp][null_number] = {"geo_potential": [], "potdrop": []}
                    
                # Initialize the dictionary for this null_number if it doesn't exist
                if null_tag not in results[timestamp][null_number]:
                    results[timestamp][null_number][null_tag] = {"geo_potential": [], "cpcp": [], "potdrop": []}
                    #results[timestamp][null_number] = {"geo_potential": [], "potdrop": []}
                
                # if null_number_plot < 80:
                #     plot_tracers_3D(ray_pls, ray_neg, nulls_pls, nulls_neg, null_number_plot)
                #     plot_tracers_2D(ray_pls, ray_neg, nulls_pls, nulls_neg, null_number_plot)
                
                # Get footprints for negative and positive traces
                Xp_sorted, Yp_sorted, Zp_sorted = get_footprints.sort_trace(ray_pls)
                Xn_sorted, Yn_sorted, Zn_sorted = get_footprints.sort_trace(ray_neg)
                
                coord_pls = get_footprints.get_footprints(Xp_sorted, Yp_sorted, Zp_sorted)
                coord_neg = get_footprints.get_footprints(Xn_sorted, Yn_sorted, Zn_sorted)
        
                # # Calculate potential drop for the positive/negative pair
                pots = get_iono_drop.get_iono_drop(iono, coord_pls, coord_neg)
    
                # # Extract the values of cpcp and potdrop from the pots dictionary
                cpcp = pots["cpcp"]
                potdrop = pots["potdrop"]
    
                # # Append the values to the lists
                results[timestamp][null_number][null_tag]["geo_potential"].append(geo_potential)
                results[timestamp][null_number][null_tag]["cpcp"].append(cpcp)
                results[timestamp][null_number][null_tag]["potdrop"].append(potdrop)
                if null_number_plot < 4:
                    plot_paneled_geo_potential(ray_neg, ray_pls, iono, nulls_neg, nulls_pls, null_number_plot, timestamp)
                #mhd_file = '../dan-files/z=0_mhd_2_e20150321-020000-000.out'                
                #plot_current(mhd_file, nulls_neg, nulls_pls)
                
                #MHD files for each timestamp for psphere runs
                mhd_file = glob.glob(f'{directory_ie}GM/z=0_mhd_2_t00{timestamp}_n*.out')
                plot_current(mhd_file[0], nulls_neg, nulls_pls, timestamp)

            except FileNotFoundError:
                print(f"Null pair not found for {filename}. Skipping...")
    
    plot_time_dynamics(results)
    # Example call to the plot function
    #mhd_file = '../dan-files/z=0_mhd_2_t00160000_n03030603.out'
    #plot_current(mhd_file, nulls_neg, nulls_pls)
    plot_mean_cpcp_vs_geo_potential(results)