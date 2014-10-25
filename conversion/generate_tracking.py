#!/usr/bin/env python
"""Save core, condensed and plume masks, 3D winds to netCDF file.

Definitions:
core - upward moving, buoyant, condensed
condensed - condensed liquid water
plume - defined following Couvreux et al. 2010
3d wind speeds - interpolated to account for Arakawa C grid

netCDF file name:
cloudtracker_input_[time].nc
"""
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset
import glob
from lib.thermo import SAM
import lib.model_param as mc

def main(filename):
    time_step = mc.time_picker(filename)
 
    # Load all the data needed to calculate core, condensed, plume and 
    # cold pools at the current time_step
    nc_file = Dataset(filename)

    tabs_field = nc_file.variables['TABS'][0,:].astype(np.double)
    qv_field = nc_file.variables['QV'][0,:].astype(np.double)/1000.
    qn_field = nc_file.variables['QN'][0,:].astype(np.double)/1000.
    p_field = nc_file.variables['p'][:].astype(np.double)*100.

    # Condensed field mask
    cloud_field = qn_field > 0.

    # Core field mask
    thetav_field = SAM.theta_v(
        p_field[:, np.newaxis, np.newaxis], tabs_field, qv_field, qn_field, 0.)
                               
    buoy_field = (thetav_field > 
         (thetav_field.mean(2).mean(1))[:, np.newaxis, np.newaxis])

    u_field = nc_file.variables['U'][0,:].astype(np.double)
    # Account for Arakawa C grid
    u_field[:, :-1, :] += u_field[:, 1:, :]
    u_field[:, -1, :] += u_field[:, 0, :]
    u_field = u_field/2.

    v_field = nc_file.variables['V'][0,:].astype(np.double)
    # Account for Arakawa C grid
    v_field[:, :, :-1] += v_field[:, :, 1:]
    v_field[:, :, -1] += v_field[:, :, 0]
    v_field = v_field/2.

    w_field = nc_file.variables['W'][0,:].astype(np.double)
    # Account for Arakawa C grid
    w_field[:-1, :, :] += w_field[1:, :, :]
    w_field[:-1, :, :] = w_field[:-1, :, :]/2.

    up_field = w_field > 0.
    core_field = up_field & buoy_field & cloud_field

    # Plume field mask; defined following Couvreux et al. 2010
    tr_field = nc_file.variables['TR01'][0, :].astype(np.double)
    x = nc_file.variables['x'][:].astype(np.double)
    y = nc_file.variables['y'][:].astype(np.double)
    z = nc_file.variables['z'][:].astype(np.double)

    nc_file.close()

    tr_mean = tr_field.reshape((len(z), len(y)*len(x))).mean(1)
    tr_stdev = np.sqrt(tr_field.reshape((len(z), len(y)*len(x))).var(1))
    tr_min = .05*np.cumsum(tr_stdev)/(np.arange(len(tr_stdev))+1)
    
    # Toggle "upward moving" constraint
    # plume_field = (tr_field > np.max(np.array([tr_mean + tr_stdev, tr_min]), \
    #    0)[:, np.newaxis, np.newaxis]) & up_field
    plume_field = (tr_field > np.max(np.array([tr_mean + tr_stdev, tr_min]), \
        0)[:, np.newaxis, np.newaxis])

    # Create netCDF file, dimensions and variables
    save_file = Dataset('%s/tracking/cloudtracker_input_%08g.nc' % \
        (mc.data_directory, time_step), 'w')
    
    save_file.createDimension('x', len(x))
    save_file.createDimension('y', len(y))
    save_file.createDimension('z', len(z))

    xvar = save_file.createVariable('x', 'f', ('x',))
    yvar = save_file.createVariable('y', 'f', ('y',))
    zvar = save_file.createVariable('z', 'f', ('z',))

    corevar = save_file.createVariable('core', 'i', ('z', 'y', 'x'))
    condvar = save_file.createVariable('condensed', 'i', ('z', 'y', 'x'))
    plumevar = save_file.createVariable('plume', 'i', ('z', 'y', 'x'))
    uvar = save_file.createVariable('u', 'f', ('z', 'y', 'x'))
    vvar = save_file.createVariable('v', 'f', ('z', 'y', 'x'))
    wvar = save_file.createVariable('w', 'f', ('z', 'y', 'x'))

    # Save dimensions and variables
    xvar[:] = x[:]
    yvar[:] = y[:]
    zvar[:] = z[:]

    corevar[:] = core_field[:]
    condvar[:] = cloud_field[:]
    plumevar[:] = plume_field[:]
    uvar[:] = u_field[:]
    vvar[:] = v_field[:]
    wvar[:] = w_field[:]
    
    save_file.close()
 
if __name__ == "__main__":
    filelist = glob.glob('%s/variables/*.nc' % mc.data_directory) 
    filelist.sort()
    nt = len(filelist)
    
    for time_step, filename in enumerate(filelist):
        print('time_step: ' + str(time_step))
        main(filename)