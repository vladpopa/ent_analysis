#!/usr/bin/env python
from pylab import *
import numpy as np
import cPickle
import glob
from netCDF4 import Dataset
from ent_analysis.lib.thermo import SAM
import var_calcs
import ent_analysis.lib.model_param as mc
import os

def index_to_zyx(index):
    ny, nx = mc.ny, mc.nx
    z = index / (ny*nx)
    index = index % (ny*nx)
    y = index / nx
    x = index % nx
    return np.array((z, y, x))

def create_savefile(t, data, vars, profile_name):
    """ Create netCDF file for specified profile at the current time step.
    netCDF axes are id and z.
    
    Return: netCDF dataset and variables
    """
    ids = data['ids'][:]
    z = data['z'][:]
    print 'cdf/%s_profile_%08d.nc' % (profile_name, t)
    savefile = Dataset(
        'cdf/%s_profile_%08d.nc' % (profile_name, t), 'w', format='NETCDF4')
    
    # Create savefile
    savefile.createDimension('ids', len(ids))
    savefile.createDimension('z', len(z))

    tsavevar = savefile.createVariable('ids', 'd', ('ids',))
    tsavevar[:] = ids[:]
    zsavevar = savefile.createVariable('z', 'd', ('z',))
    zsavevar[:] = z[:]

    variables = {}
    for name in vars:
        variables[name] = savefile.createVariable(name, 'd', ('ids', 'z'))
        
    return savefile, variables

def make_profile(z_indexes, y_indexes, x_indexes, data, vars, profiles):
    z = np.unique(z_indexes)
    # Average properties at each height
    for k in z:
        mask = (z_indexes == k)
        if mask.sum() == 0: continue
        j = y_indexes[mask]
        i = x_indexes[mask]
        for name in vars:
            profiles[name][k] = vars[name](data, k, j, i)

    return profiles

def make_profiles(profiles, cloud_data, vars, data, n):
    """Make profiles for a single cloud at the current time step.
    """
    # Iterate through profile types for a single cloud
    for item in ('core', 'condensed', 'plume',
        'condensed_shell', 'condensed_edge', 'condensed_env',
        'core_shell', 'core_edge', 'core_env', 
        'plume_shell', 'plume_edge', 'plume_env'):  
           
        variables = profiles[item]                
        temp_profile = {}
        # Initialize empty profile for all variables for a profile type
        for name in vars:
            temp_profile[name] = ones_like(data['z'][:])*np.NaN

        # Make profile if it contains any mask points 
        indexes = cloud_data[item]
        if len(indexes) > 0:
            z, y, x = index_to_zyx(indexes)            
            results = make_profile(z, y, x, data, vars, temp_profile)
        else:
            results = temp_profile       
              
        # Index by profile type and index (not id) of current cloud                                 
        for name in vars:
            variables[name][n, :] = results[name]

def main(filename):
    """Create profiles at current time step.
    """
    vars = {
          'AREA': var_calcs.area,
          'TABS': var_calcs.tabs,
          'QN': var_calcs.qn,
          'QV': var_calcs.qv,
          'QT': var_calcs.qt,
          'U': var_calcs.u,
          'V': var_calcs.v,
          'W': var_calcs.w,
          'THETAV': var_calcs.thetav,
          'THETAV_LAPSE': var_calcs.thetav_lapse,
          'THETAL': var_calcs.thetal,
          'MSE': var_calcs.mse,
          'RHO': var_calcs.rho,
          'PRES': var_calcs.press,
          'WQREYN': var_calcs.wqreyn,
          'WWREYN': var_calcs.wwreyn,
          'DWDZ': var_calcs.dw_dz,
          'DPDZ': var_calcs.dp_dz,
          'TR01': var_calcs.tr01,
          'RELH': var_calcs.relh
    }
    
    # Automatically load time step from output file name
    time = mc.time_picker(filename)
    
    # Load netCDF files
    nc_file = Dataset(filename)
    stat_file = Dataset(mc.get_stat())

    data = {'z': nc_file.variables['z'][:].astype(double),
            'p': nc_file.variables['p'][:].astype(double),
            'RHO' : stat_file.variables['RHO'][time,:].astype(double),
            }
    stat_file.close()
   
    # Load the cloud data at current timestep
    cloud_filename = '../cloudtracker/pkl/cloud_data_%08d.pkl' % time
    clouds = cPickle.load(open(cloud_filename, 'rb'))
       
    ids = clouds.keys()
    ids.sort()

    data['ids'] = np.array(ids)
    for name in ('QV', 'QN', 'TABS', 'PP', 'U', 'V', 'W', 'TR01'):
        data[name] = nc_file.variables[name][0, :].astype(np.double)
                
    # Create a savefile for each profile at the current time step
    savefiles = {}
    profiles = {}
    for item in ('core', 'condensed', 'plume',
        'condensed_shell', 'condensed_edge', 'condensed_env',
        'core_shell', 'core_edge', 'core_env', 
        'plume_shell', 'plume_edge', 'plume_env'):            
		 
        savefile, variables = create_savefile(time, data, vars, item)
        savefiles[item] = savefile
        profiles[item] = variables
        
    for n, id in enumerate(ids):
        print "time: ", time, " id: ", id
        # Select the current cloud id
        cloud = clouds[id]
	
        make_profiles(profiles, cloud, vars, data, n)
        
    for savefile in savefiles.values():
        savefile.close()
    nc_file.close()
   
if __name__ == "__main__":
    files = glob.glob('%s/variables/*.nc' % mc.data_directory)
    files.sort()
    
    for time, filename in enumerate(files):
        main(filename)