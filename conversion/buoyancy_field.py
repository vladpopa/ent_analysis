#!/usr/bin/env python
"""Save buoyancy field to netCDF file.

netCDF file name:
buoyancy_[time].nc
"""
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset
import glob
from lib.thermo import SAM
import lib.model_param as mc

def calculate_buoyancy(filename):
    time_step = mc.time_picker(filename)
 
    # Load all the data needed to calculate buoyancy field at current time step
    nc_file = Dataset(filename)

    tabs_field = nc_file.variables['TABS'][0,:].astype(np.double)
    qv_field = nc_file.variables['QV'][0,:].astype(np.double)/1000.
    qn_field = nc_file.variables['QN'][0,:].astype(np.double)/1000.
    p_field = nc_file.variables['p'][:].astype(np.double)*100.

    # Core field mask
    thetav_field = SAM.theta_v(
        p_field[:, np.newaxis, np.newaxis], tabs_field, qv_field, qn_field, 0.)

    x = nc_file.variables['x'][:].astype(np.double)
    y = nc_file.variables['y'][:].astype(np.double)
    z = nc_file.variables['z'][:].astype(np.double)

    nc_file.close()

    # Calculate buoyanct field
    buoy_field = SAM.g*(thetav_field - \
        (thetav_field.mean(2).mean(1))[:, np.newaxis, np.newaxis])/ \
        (thetav_field.mean(2).mean(1))[:, np.newaxis, np.newaxis]

    # Create netCDF file, dimensions and variables
    save_file = Dataset(
        '%s/tracking/buoyancy_%08g.nc' % (mc.data_directory, time_step), 'w')
    
    save_file.createDimension('x', len(x))
    save_file.createDimension('y', len(y))
    save_file.createDimension('z', len(z))

    xvar = save_file.createVariable('x', 'f', ('x',))
    yvar = save_file.createVariable('y', 'f', ('y',))
    zvar = save_file.createVariable('z', 'f', ('z',))
    
    xvar.units = 'm'
    xvar.long_name = 'x coordinate'
    yvar.units = 'm'
    yvar.long_name = 'y coordinate'
    zvar.units = 'm'
    zvar.long_name = 'z coordinate'
        
    buoyvar = save_file.createVariable('B', 'f', ('z', 'y', 'x'))
    buoyvar.units = 'm/s2'
    buoyvar.long_name = 'Buoyancy'
    
    # Save dimensions and variables
    xvar[:] = x[:]
    yvar[:] = y[:]
    zvar[:] = z[:]

    buoyvar[:] = buoy_field[:]
    
    save_file.close()
 
if __name__ == "__main__":
    filelist = glob.glob('%s/variables/*.nc' % mc.data_directory) 
    filelist.sort()
    nt = len(filelist)
    
    for time_step, filename in enumerate(filelist):
        print('time_step: ' + str(time_step))
        calculate_buoyancy(filename)