#!/usr/bin/env python

"""
Calculate cloud lifetime, mean cloud mass, maximum cloud height and minimum
cloud base for tracked clouds.

Examples
-------
$ python cloud_statistics.py
"""

from __future__ import division, print_function
import sys
sys.path.append('../lib')
import glob
import os
import numpy as np
import numpy.ma as ma
import model_param as mc
from netCDF4 import Dataset
   
def cloud_statistics(file_name):
    """
    Return cloud lifetime, mean cloud mass, maximum cloud height and minimum
    cloud base for a tracked cloud.
        
    Parameters
    ----------
    file_name : netCDF file name
        id_profile file for a tracked cloud with dimensions double t(t), 
        double z(z).
      
    Return
    ------
    numpy.array
        Column array of shaded cloud areas.    
    """
    
    # Read netCDF dataset
    data = Dataset(file_name)
    
    # Calculate cloud lifetime in seconds
    time = data.variables['t'][...]
    lifetime = len(time)*mc.dt

    # Calculate cloud base and cloud top in metres
    area = ma.masked_invalid(data.variables['AREA'][...])
    z = ma.masked_array(data.variables['z'][...], np.all(ma.getmask(area), \
        axis=0))
    base = z.min()
    top = z.max()

    # Calculate cloud mass in kilograms
    qn = ma.masked_invalid(data.variables['QN'][...])
    rho = ma.masked_invalid(data.variables['RHO'][...])
    mass = np.mean(np.sum(area*rho*mc.dz, axis=1))
    # Need to add QV and QN here
    
    return lifetime, base, top, mass

if __name__ == '__main__':
      
    # Create pkl directory to store cloud statistics
    if not os.path.exists('npy'):
        os.makedirs('npy')
    
    # Collect all netCDF files containing cloud id_profiles
    input_files = os.path.join('../id_profiles/cdf', 'condensed_profile_*.nc')
    file_list = glob.glob(input_files)
    file_list.sort()
    
    # Calculate and store statistics for all tracked clouds
    # stats = ['max_height', 'min_height', 'average_mass', 'duration']
    # cloud_stat = dict((key, np.array([])) for key in stats)
    
    lifetimes = np.array([])
    bases = np.array([])
    tops = np.array([])    
    masses = np.array([])   
    for n, file_name in enumerate(file_list):
        print("Calculating statistics for cloud id %d" % \
             int(file_name[-11:-3]))
        lifetime, base, top, mass = cloud_statistics(file_name)     
        lifetimes = np.append(lifetimes, np.array(lifetime))
        bases = np.append(bases, base)
        tops = np.append(tops, top)
        masses = np.append(masses, mass)
        
    # Discard cloud noise?
    # cloud_stat['duration'] =  lifetimes[1:]
    # cloud_stat['min_height'] =  bases[1:]
    # cloud_stat['max_height'] =  tops[1:]
    
    # Save cloud statistics
    np.save('npy/cloud_statistics_lifetimes.npy', lifetimes[1:])
    np.save('npy/cloud_statistics_bases.npy', bases[1:])
    np.save('npy/cloud_statistics_tops.npy', tops[1:])
    np.save('npy/cloud_statistics_masses.npy', masses[1:])