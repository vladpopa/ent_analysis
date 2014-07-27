#!/usr/bin/env python

"""
Calculate core duration, minimum core base, maximum core height, mean core mass,
formation time, dissipation time, maximum depth, depth evolution and 
corresponding times for tracked clouds.

Output is saved in a pkl file as a list of dictionaries, with one dictionary per
core id. Dictionary keys are: 'id', 'duration', 'min_height', 'max_height', 
'average_mass', 'l_min', 'l_max', 'depth', 'max_depth', 'time'].

Examples
-------
$ python core_id_stats.py
"""

from __future__ import division, print_function
import glob
import os
import numpy as np
import numpy.ma as ma
import model_config.cgils_ctl_s6_25m as mc
from netCDF4 import Dataset
try:
    import cPickle as pickle
except:
    import pickle     
       
def cloud_statistics(file_name):
    """
    Return core duration, minimum core base, maximum core height, mean core 
    mass, formation time, dissipation time, maximum depth, depth evolution and 
    corresponding times for tracked clouds.
        
    Parameters
    ----------
    file_name : netCDF file name
        id_profile file for a tracked core with dimensions double t(t), 
        double z(z).
      
    Return
    ------
    tuple : cloud_id, lifetime, base, top, mass, l_min, l_max, depths,
        max_depth, times
    """
    
    # Read netCDF dataset
    data = Dataset(file_name)
    
    # Core ID
    cloud_id = int(file_name[-11:-3])
    
    # Core duration (seconds)
    times = data.variables['t'][...]
    lifetime = len(times)*mc.dt
    
    # Formation time, dissipation time (seconds)
    l_min = times.min()*mc.dt
    l_max = times.max()*mc.dt

    # Minimum core base, maximum core height, maximum depth, depth evolution 
    # (metres)
    area = ma.masked_invalid(data.variables['AREA'][...])
    z = data.variables['z'][...]
    z = z*np.ones(np.shape(area))
    z = ma.masked_array(z, ma.getmask(area)) 
    bases = z.min(axis=1)
    tops = z.max(axis=1)
    depths = tops - bases + mc.dz
    max_depth = depths.max()
    base = bases.min()
    top = tops.max()

    # Mean core mass mass (kilograms)
    qn = ma.masked_invalid(data.variables['QN'][...])
    rho = ma.masked_invalid(data.variables['RHO'][...])
    mass = np.mean(np.sum(area*rho*mc.dz, axis=1))

    # Remove missing values
    times = ma.masked_array(times, ma.getmask(depths))
    depths = depths[~depths.mask]
    times = times[~times.mask]
    
    data.close()
    
    return cloud_id, lifetime, base, top, mass, l_min, l_max, depths, \
        max_depth, times

if __name__ == '__main__':    
    # Create pkl directory to store core statistics
    if not os.path.exists('pkl'):
        os.makedirs('pkl')
        
    # Collect all netCDF files containing core profiles
    input_files = os.path.join('../id_profiles/cdf', 'core_profile_*.nc')
    file_list = glob.glob(input_files)
    file_list.sort()
    
    # Calculate and store core statistics for all tracked clouds
    keys = ['id', 'duration', 'min_height', 'max_height', 'average_mass', 
        'l_min', 'l_max', 'depth', 'max_depth', 'time']
    id_stats = []
    
    for n, file_name in enumerate(file_list):
        print("Calculating statistics for cloud id %d" % int(file_name[-11:-3]))
        id_stats.append(dict(zip(keys, cloud_statistics(file_name))))
        
    # Find id of longest lived core
    durations = np.array([])
    ids = np.array([])
    for item in id_stats:
        durations = np.hstack((durations, item['depth']))
        ids = np.hstack((ids, item['id']))
    idx_max = np.argmax(durations)
    print('max_lifetime:', int(durations[idx_max]), 'id:', int(ids[idx_max]))
    
    # Save core statistics
    pickle.dump(id_stats, open('pkl/core_id_stats.pkl', 'wb'))