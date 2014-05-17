#!/usr/bin/env python

"""
Calculate cloud lifetime, mean cloud mass, maximum cloud height and minimum
cloud base for tracked clouds.

Examples
-------
$ python cloud_statistics.py
"""

from __future__ import division, print_function
import glob
import os
import numpy as np
import numpy.ma as ma
import ent_analysis.lib.model_param as mc
from netCDF4 import Dataset
try:
    import cPickle as pickle
except:
    import pickle
       
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
    l_min = time.min()*mc.dt
    l_max = time.max()*mc.dt

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
    
    data.close()
    
    return lifetime, base, top, mass, l_min, l_max

if __name__ == '__main__':
      
    # Create pkl directory to store cloud statistics
    if not os.path.exists('pkl'):
        os.makedirs('pkl')
    
    # Collect all netCDF files containing cloud id_profiles
    input_files = os.path.join('../id_profiles/cdf', 'condensed_profile_*.nc')
    file_list = glob.glob(input_files)
    file_list.sort()
    
    # Calculate and store statistics for all tracked clouds
    stats = ['max_height', 'min_height', 'average_mass', 'duration', 'l_min' \
        'l_max']
    id_stats = dict((key, np.array([])) for key in stats)
    
    lifetimes = np.array([])
    bases = np.array([])
    tops = np.array([])    
    masses = np.array([]) 
    l_mins = np.array([]) 
    l_maxs = np.array([])    
    ids = np.array([])     
    for n, file_name in enumerate(file_list):
        print("Calculating statistics for cloud id %d" % \
             int(file_name[-11:-3]))
        lifetime, base, top, mass, l_min, l_max = cloud_statistics(file_name)     
        lifetimes = np.append(lifetimes, np.array(lifetime))
        bases = np.append(bases, base)
        tops = np.append(tops, top)
        masses = np.append(masses, mass)
        l_mins = np.append(l_mins, l_min)
        l_maxs = np.append(l_maxs, l_max)
        ids = np.append(ids, int(file_name[-11:-3]))
        
    # Remove cloud noise
    id_stats['duration'] =  lifetimes[1:]/60.
    id_stats['min_height'] =  bases[1:]
    id_stats['max_height'] =  tops[1:]
    id_stats['l_min'] = l_mins[1:]/60.
    id_stats['l_max'] = l_maxs[1:]/60.
    id_stats['average_mass'] = masses[1:]
    id_stats['id'] = ids[1:]
    
    # Find id of longest lived cloud
    id_max = id_stats['id'][np.argmax(id_stats['duration'])]
    print('max_lifetime:', id_stats['duration'][id_max], 'id:', id_max)
    
    # Save cloud statistics
    pickle.dump(id_stats, open('pkl/stats.pkl', 'wb'))