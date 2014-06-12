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
import bomex  as mc
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
    # Need to add QV and QN here ??? or just dry air mass???

    data.close()
    
    return lifetime, base, top, mass, l_min, l_max

def subcloudplume_statistics(file_name, base):
    """
    Return sub-cloud plume lifetime for a tracked cloud.
        
    Parameters
    ----------
    file_name : netCDF file name
        id_profile file for a tracked cloud with dimensions double t(t), 
        double z(z).
    base : cloud base (in m)
      
    """
    # Read netCDF dataset
    data = Dataset(file_name)
    tracer = ma.masked_invalid(data.variables['TR01'][...])
    area = ma.masked_invalid(data.variables['AREA'][...])
    time = data.variables['t'][...]
    z = data.variables['z'][...]
    
    # find grid points above cloud base
    above_cloudbase = z >=  base
    above_cloudbase = np.tile(above_cloudbase, (time.size, 1))

    # set this to define the lifetime of a sub-cloud plume
    # if the maximum tracer concentration below cloudbase is equal to
    # or less than this value, then the plume has dissipated
    threshold = 1

    # mask is True for points that don't have tracer
    tracer_mask = ma.getmask(tracer)
    
    area_mask = ma.getmask(area)

    tracer_below_cloudbase = ma.masked_array(tracer, np.logical_or(above_cloudbase, tracer_mask))
    tracer_max = np.max(tracer_below_cloudbase, axis = 1)

    area_below_cloudbase = ma.masked_array(area, np.logical_or(above_cloudbase, area_mask))
    area_mean = np.mean(area_below_cloudbase, axis = 1)
    area_mean = np.mean(area_mean)
 
    # if there is no tracer below cloudbase, the subcloud plume has zero lifetime
    if np.all(ma.getmask(tracer_max)):
       lifetime  = 0
       l_min = np.min(time)
       l_max = np.min(time)
    
    #otherwise get the time when the max. tracer concentration below cloudbase falls below the threshold
    else:
        tend = time[tracer_max <= threshold]
        if len(tend):
            tend = np.min(tend)
            tstart = np.min(time)
            l_min = tstart*mc.dt
            l_max = tend*mc.dt 
            lifetime = (tend-tstart)*mc.dt
        else:
            lifetime = len(time)*mc.dt
            l_min = np.min(time)*mc.dt
            l_max = np.max(time)*mc.dt


    data.close()

    return lifetime, l_min, l_max, area_mean
    
if __name__ == '__main__':
      
    # Create pkl directory to store cloud statistics
    if not os.path.exists('pkl'):
        os.makedirs('pkl')
    
    # Collect all netCDF files containing cloud id_profiles
    plume_input_files = os.path.join(mc.analysis_dir, 'id_profiles/cdf', \
        'plume_profile_*.nc')
    plume_file_list = glob.glob(plume_input_files)
    plume_file_list.sort()

    
    # Calculate and store statistics for all tracked clouds
    stats = ['max_height', 'min_height', 'average_mass', 'cloud_duration', 'l_min' \
        'l_max']
    id_stats = dict((key, np.array([])) for key in stats)
    
    cloud_lifetimes = np.array([])
    bases = np.array([])
    tops = np.array([])    
    masses = np.array([]) 
    l_mins = np.array([]) 
    l_maxs = np.array([])    
    cloud_ids = np.array([])
    plume_lifetimes = np.array([])
    plume_l_mins = np.array([])
    plume_l_maxs = np.array([])
    plume_areas = np.array([])
    
    for n, plume_file_name in enumerate(plume_file_list):
        id = int(plume_file_name[-11:-3])
        cloud_file_name = glob.glob('id_profiles/cdf/condensed_profile_%08d.nc' % id)
        print('Calculating statistics for plume id', id)
        cloud_lifetime, base, top, mass, l_min, l_max = cloud_statistics(cloud_file_name[0]) 
        plume_lifetime, plume_l_min, plume_l_max, plume_area = subcloudplume_statistics(plume_file_name, base)
        cloud_lifetimes = np.append(cloud_lifetimes, cloud_lifetime)
        bases = np.append(bases, base)
        tops = np.append(tops, top)
        masses = np.append(masses, mass)
        l_mins = np.append(l_mins, l_min)
        l_maxs = np.append(l_maxs, l_max)
        plume_lifetimes = np.append(plume_lifetimes, plume_lifetime) 
        plume_l_mins = np.append(plume_l_mins, plume_l_min)
        plume_l_maxs = np.append(plume_l_maxs, plume_l_max)
        plume_areas = np.append(plume_areas, plume_area)
        cloud_ids = np.append(cloud_ids, id)
        print('cloud duration', cloud_lifetime/mc.dt)
        print('sub-cloud plume duration', plume_lifetime/mc.dt)
        print('sub-cloud plume mean area', plume_area)
        
    # Remove cloud noise
    id_stats['cloud_duration'] =  cloud_lifetimes[1:]/mc.dt
    id_stats['min_height'] =  bases[1:]
    id_stats['max_height'] =  tops[1:]
    id_stats['l_min'] = l_mins[1:]/mc.dt
    id_stats['l_max'] = l_maxs[1:]/mc.dt
    id_stats['average_mass'] = masses[1:]
    id_stats['cloud_id'] = cloud_ids[1:]
    id_stats['plume_duration'] = plume_lifetimes[1:]/mc.dt
    id_stats['plume_l_min'] = plume_l_mins[1:]/mc.dt
    id_stats['plume_l_max'] = plume_l_maxs[1:]/mc.dt
    id_stats['plume_meanarea'] = plume_areas[1:]
 

         
    # Find id of longest lived cloud
    id_max = id_stats['cloud_id'][np.argmax(id_stats['cloud_duration'])]
    print('max_lifetime:', id_stats['cloud_duration'][id_max], 'id:', id_max)
    print('mean area of sub-cloud plume  of longest lived cloud', id_stats['plume_meanarea'][id_max])

    print('mean cloud lifetime (sec): ', np.mean(id_stats['cloud_duration'])*mc.dt)
    print('mean subcloud plume lifetime (sec): ', np.mean(id_stats['plume_duration'])*mc.dt)
  

    
    # Save cloud statistics
    pickle.dump(id_stats, open('pkl/stats_test.pkl', 'wb'))

