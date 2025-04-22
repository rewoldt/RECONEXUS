#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:07:53 2024

@author: rewoldt
"""

from glob import glob
import numpy as np
    
import os
import re
from datetime import datetime, timedelta

def read_start_time(param_file_path):
    '''
    Reads start time from PARAM.in fileand returns the date time in %Y%m%d-%H%M%S.
    
    Parameters
    ----------
    
    '''
    with open(param_file_path, 'r') as f:
        lines = f.readlines()
    
    # Find the start time block
    for i, line in enumerate(lines):
        if "#STARTTIME" in line:
            year = int(lines[i+1].strip()[0:4])
            month = int(lines[i+2].strip()[1])
            day = int(lines[i+3].strip()[1])
            hour = int(lines[i+4].strip()[0])
            minute = int(lines[i+5].strip()[0])
            second = int(lines[i+6].strip()[0])
            return datetime(year, month, day, hour, minute, second)
    
    raise ValueError("Start time not found in PARAM.in file")

def extract_time_from_filename(filename):
    '''
    Extracts the time from the filename in the format 'tHHMMSS'.
    '''
    match = re.search(r'_t(\d{8})\.dat$', filename)
    if match:
        time_str = match.group(1)
        print (time_str)
        hours = int(time_str[0:4])
        print (hours)
        minutes = int(time_str[4:6])
        print (minutes)
        seconds = int(time_str[6:8])
        print (seconds)
        return timedelta(hours=hours, minutes=minutes, seconds=seconds)
    return None

def rename_files_in_directory(directory, start_time):
    '''
    Renames all files matching the patterns by updating their time.
    '''
    for filename in os.listdir(directory):
        # Check if the file matches the pattern
        if re.match(r'(null_line_neg.*_t\d{8}\.dat|null_line_pls.*_t\d{8}\.dat|NegNulls_t\d{8}\.dat|PlusNulls_t\d{8}\.dat)', filename):
            # Extract the time from the filename
            time_delta = extract_time_from_filename(filename)
            
            #if time_delta is not None:
            # Calculate the new time by adding the delta to the start time
            new_time = start_time + time_delta
            new_time_str = new_time.strftime('%Y%m%d-%H%M%S')
            
            # Form the new filename
            new_filename = re.sub(r'_t\d{8}', f'_e{new_time_str}', filename)
            
            # Rename the file
            old_path = os.path.join(directory, filename)
            new_path = os.path.join(directory, new_filename)
            os.rename(old_path, new_path)
            print(f'Renamed {filename} to {new_filename}')

if __name__ == "__main__":
    # Directory containing the files
    param_directory = "/Volumes/coe-dwelling/swmf_results/GmPlasSphere/plas_psgmie_smallby"
    run_directory = "/Volumes/SWMF_runs/reconnection_perfection/RECONX/run"
    #run_directory = "/Volumes/coe-dwelling/swmf_results/reconnection_perfection/RECONX/run"
    
    # Path to PARAM.in
    param_file = os.path.join(param_directory, "PARAM.in.ss")
    
    # Read the start time
    start_time = read_start_time(param_file)
    
    # Rename the files
    rename_files_in_directory(run_directory, start_time)
