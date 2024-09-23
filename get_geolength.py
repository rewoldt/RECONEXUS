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
from spacepy.pybats import rim
from datetime import datetime
from scipy.stats import pearsonr


def cross_product(u, v):
    dim = len(u)
    s = []
    for i in range(dim):
        if i == 0:
            j, k = 1, 2
            s.append(u[j]*v[k] - u[k]*v[j])
        elif i == 1:
            j, k = 2, 0
            s.append(u[j]*v[k] - u[k]*v[j])
        else:
            j, k = 0, 1
            s.append(u[j]*v[k] - u[k]*v[j])
    return s


def dot_product(u, v):
    return sum([un*vn for un, vn in zip(u, v)])


def calc_uxb(u, B):
    uxb = cross_product(u, B)
    return uxb


def get_potential(null_file_neg, null_file_pls, n):

    Xp_sorted, Yp_sorted, Zp_sorted = get_footprints.sort_trace(null_file_pls)
    Xn_sorted, Yn_sorted, Zn_sorted = get_footprints.sort_trace(null_file_neg)
    # print (Xn_sorted, Yn_sorted, Zn_sorted)

    Xn = Xn_sorted[n:]
    Yn = Yn_sorted[n:]
    Zn = Zn_sorted[n:]
    Xp = Xp_sorted[n:]
    Yp = Yp_sorted[n:]
    Zp = Zp_sorted[n:]

    RE2km = 6378.14  # *(10**3) #km

    # geo_potential = []
    dist = []
    length = []
    length_dict = {}

    uxb = calc_uxb(u, B)

    for i in range(len(Xn)):
        for j in range(len(Xp)):
            distance = (np.sqrt((Xn[i]-Xp[j])**2 +
                                (Yn[i]-Yp[j])**2 +
                                (Zn[i]-Zp[j])**2)*RE2km)
            dist.append(distance)
            L = (Xn[i]-Xp[j])*RE2km, (Yn[i]-Yp[j])*RE2km, (Zn[i]-Zp[j])*RE2km
            length.append([(Xn[i]-Xp[j])*RE2km,
                           (Yn[i]-Yp[j])*RE2km,
                           (Zn[i]-Zp[j])*RE2km])

    min_distance_index = np.argmin(dist)
    min_distance = dist[min_distance_index]

    #print (min_distance_index, min_distance, length[min_distance_index],
    #       length[min_distance_index][0], length[min_distance_index][1],
    #       length[min_distance_index][2],
    #       (np.sqrt((length[min_distance_index][0])**2 +
    #       (length[min_distance_index][1])**2 +
    #       (length[min_distance_index][2])**2)))

    geo_potential = dot_product(uxb, length[min_distance_index])

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
                    times.append(datetime.strptime(timestamp, '%H%M%S'))
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
    plt.scatter(mean_geo_potential_values, mean_cpcp_values, label=f"Null Number {null_number}", marker='^', color='b')

        # Calculate correlation coefficient
    correlation, _ = pearsonr(mean_geo_potential_values, mean_cpcp_values)

    # Add text for correlation
    plt.text(0.05, 0.95, f'Correlation: {correlation:.2f}', ha='left', va='center', transform=plt.gca().transAxes)
    plt.xscale('log')
    plt.yscale('log')
    #plt.semilogy(times, mean_cpcp_values, label=f'CPCP', marker='^', linestyle=':', color='b')
    #plt.semilogy(times, mean_potdrop_values, label=f'Potential drop', marker='s', linestyle='--', color='g')
    #plt.semilogy(times, mean_geo_potential_values, label=f'GEL Potential', marker='o', color='r')
    # Customize the plot
    plt.xlabel('GEL Potential')
    plt.ylabel('CPCP (kV)')
    plt.title('Comparison of  CPCP and GEL Potential for All Null Pairs')
    plt.legend(loc='best', ncol=3)#, prop={'size': 5.5})
    plt.grid(True)
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save the plot if needed
    plt.savefig('all_nulls_mean_comparison_plot.png', dpi=300)

    # Show the plot
    plt.show()

def plot_mean_cpcp_vs_geo_potential(results):
    mean_cpcp_values = []
    mean_geo_potential_values = []
    mean_potdrop_values = []
    null_pairs = []

    for timestamp in results:
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

    # Plotting Mean CPCP vs Mean Geo Potential
    #plt.scatter(mean_geo_potential_values, mean_cpcp_values, marker='8', color='b', s=80, label='CPCP vs GEL Potential')
    plt.scatter(mean_geo_potential_values, mean_cpcp_values, s=70, label='CPCP vs GEL Potential')
    # Calculate correlation coefficient for Mean CPCP vs Mean Geo Potential
    #correlation_cpcp, _ = pearsonr(mean_geo_potential_values, mean_cpcp_values)
    # Add text for correlation for Mean CPCP vs Mean Geo Potential
    #plt.text(0.05, 0.95, f'CPCP Correlation: {correlation_cpcp:.2f}', ha='left', va='center', transform=plt.gca().transAxes)

    # Plotting Mean Potdrop vs Mean Geo Potential
    #plt.scatter(mean_geo_potential_values, mean_potdrop_values, marker='8', color='r', s=80, label='Iono Footprint Potential Drop vs GEL Potential')
    plt.scatter(mean_geo_potential_values, mean_potdrop_values, s=70, color='indianred', label='Iono Footprint Potential Drop vs GEL Potential')

    # Calculate correlation coefficient for Mean Potdrop vs Mean Geo Potential
    #correlation_potdrop, _ = pearsonr(mean_geo_potential_values, mean_potdrop_values)
    # Add text for correlation for Mean Potdrop vs Mean Geo Potential
    #plt.text(0.05, 0.90, f'Potdrop Correlation: {correlation_potdrop:.2f}', ha='left', va='center', transform=plt.gca().transAxes)

    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('GEL Potential (kV)', fontsize=16)
    plt.ylabel('Ionospheric Potential (kV)', fontsize=16)
    plt.title('Iono Footprint Potential Drop and CPCP vs. GEL Potential', fontsize=18)
    plt.grid(True)
    plt.legend(fontsize=14)
    plt.tight_layout()

    # Save the plot
    plt.savefig('potential_comparison.png', dpi=300)

    # Show the plot
    plt.show()

if __name__ == '__main__':
    # Your initial code goes here to populate the results dictionary
    ux = -400 * 10**3  # m/s
    uy = 0
    uz = 0

    Bx = 0
    By = -2 * 10**-9  # T
    Bz = -10 * 10**-9  # T

    # ux = -450 * 10**3  # m/s
    # uy = 0
    # uz = 0

    # Bx = 0
    # By = 2 * 10**-9  # T
    # Bz = -10 * 10**-9  # T

    B = [Bx, By, Bz]
    u = [ux, uy, uz]

    uxb = calc_uxb(u, B)

    # Open a file:
    #iono = rim.Iono("/Users/rewoldt/HighResAMRdemo-ReconxWork/IE/it150321_000100_000.idl")
    # iono = {}
    # directory_ie = "/Users/rewoldt/HighResAMRdemo-ReconxWork/IE/"
    # #directory_ie="/Users/rewoldt/ideal_recon_tim/IE/"
    # for filename in os.listdir(directory_ie):
    #     #if filename.endswith(".tec"):
    #     if filename.endswith(".idl"):
    #         # Extract timestamp from the filename
    #         timestamp = filename.split("_")[1]
    #         iono[timestamp] = rim.Iono(os.path.join(directory_ie, filename))

    # Define the directory containing the trace files
    directory = "/Users/rewoldt/RECONX/run/"
    #directory = "/Users/rewoldt/RECONX/run_amr/"
    pattern = re.compile(r"null_line_neg_o\d+_\d{3}_e\d+-\d{6}\.dat")
    #pattern = re.compile(r"null_line_neg_o\d+_\d{3}_n\d+.dat") #null_line_neg_o03_001_n03030603.dat
    #pattern = re.compile(r"null_line_neg_[ns]\d+_\d{3}_e\d+-\d{6}\.dat")
    #pattern = re.compile(r"null_line_neg_[ns]\d+_001_e\d+-\d{6}\.dat")
    results = {}


    # Loop over each trace file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".dat") and pattern.match(filename):
            # Extract timestamp from the filename
            timestamp = filename.split("_")[-1].split("-")[-1].split(".")[0]

            # Extract null number and polarity from the filename
            null_number = filename.split("_")[3][1:]
            polarity = filename.split("_")[3][:1]
            null_tag = filename.split("_")[4]
            nulls_neg = directory+f"NegNulls_e20150321-{timestamp}.dat"
            nulls_pls = directory+f"PlusNulls_e20150321-{timestamp}.dat"
            #nulls_neg = directory+f"NegNulls_n03030603.dat"
            #nulls_pls = directory+f"PlusNulls_n03030603.dat"
            null_number_plot = int(null_number) - 1
            print (null_number, null_number_plot)

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

                if null_number_plot < 80:
                    plot_tracers_3D(ray_pls, ray_neg, nulls_pls, nulls_neg, null_number_plot)
                    plot_tracers_2D(ray_pls, ray_neg, nulls_pls, nulls_neg, null_number_plot)

                # Get footprints for negative and positive traces
                # Xp_sorted, Yp_sorted, Zp_sorted = get_footprints.sort_trace(ray_pls)
                # Xn_sorted, Yn_sorted, Zn_sorted = get_footprints.sort_trace(ray_neg)

                # coord_pls = get_footprints.get_footprints(Xp_sorted, Yp_sorted, Zp_sorted)
                # coord_neg = get_footprints.get_footprints(Xn_sorted, Yn_sorted, Zn_sorted)

                # # Calculate potential drop for the positive/negative pair
                # pots = get_iono_drop.get_iono_drop(iono[timestamp], coord_pls, coord_neg)

                # # Extract the values of cpcp and potdrop from the pots dictionary
                # cpcp = pots["cpcp"]
                # potdrop = pots["potdrop"]

                # # Append the values to the lists
                results[timestamp][null_number][null_tag]["geo_potential"].append(geo_potential)
                # results[timestamp][null_number][null_tag]["cpcp"].append(cpcp)
                # results[timestamp][null_number][null_tag]["potdrop"].append(potdrop)

            except FileNotFoundError:
                print(f"Null pair not found for {filename}. Skipping...")

    # Example call to the plot function
    #plot_results(results)
    plot_mean_cpcp_vs_geo_potential(results)
