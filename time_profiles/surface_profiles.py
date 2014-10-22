#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import cPickle
from netCDF4 import Dataset
import os
from lib.thermo import SAM
import lib.model_param as mc

# Load mean cloud field statistics
stat_file = Dataset(mc.get_stat())
data = {'z': stat_file.variables['z'][:].astype(np.double),
    'RHO' : stat_file.variables['RHO'][0,:].astype(np.double),
    'PRES' : stat_file.variables['PRES'][0,:].astype(np.double)*100.}
stat_file.close()

def create_savefile(t, data, vars, profile_name):
    """ Create netCDF file for specified profile at the current time step.
    netCDF axes are id and z.
    
    Return: netCDF dataset and variables
    """
    ids = data['ids'][:]
    z = data['z'][:]
    savefile = Dataset('cdf/%s_profile_%08d.nc' % (profile_name, t), 'w')
    
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

def make_profiles(variables, cloud_data, vars, data, n):
    """Make profiles for a single cloud at the current time step.
    """
    temp = {'core': 'CORE_SURFACE', 'condensed': 'CONDENSED_SURFACE', 
        'plume': 'PLUME_SURFACE'}
    
    for item in temp:                
        temp_profile = {}
        for name in vars:
            temp_profile[name] = np.zeros_like(data['z'][:])

        indexes = cloud_data[item]
        if len(indexes) > 0:
            # For masks of top and bottom cloud surface
            indexes2 = indexes + mc.nx*mc.ny
            a = ~np.in1d(indexes, indexes2, assume_unique=True)
            indexes2 = indexes - mc.nx*mc.ny
            b = ~np.in1d(indexes, indexes2, assume_unique=True)
        
            z, y, x = mc.index_to_zyx(indexes)
        
            # Find masks of lateral boundaries in x direction; account for 
            # reentrant domain
            x2  = (x+1) % mc.nx
            indexes2 = mc.zyx_to_index(z, y, x2)
            c = ~np.in1d(indexes, indexes2, assume_unique=True)
            x2  = (x-1) % mc.nx
            indexes2 = mc.zyx_to_index(z, y, x2)
            d = ~np.in1d(indexes, indexes2, assume_unique=True)

            # Find masks of lateral boundaries in y direction; account for 
            # reentrant domain
            y2  = (y+1) % mc.ny
            indexes2 = mc.zyx_to_index(z, y2, x)
            e = ~np.in1d(indexes, indexes2, assume_unique=True)
            y2  = (y-1) % mc.ny
            indexes2 = mc.zyx_to_index(z, y2, x)
            f = ~np.in1d(indexes, indexes2, assume_unique=True)
            
            # Mask of cloud surface converted to dimensional units
            area = mc.dx*mc.dy*(a+b) + mc.dy*mc.dz*(c+d) + mc.dx*mc.dz*(e+f)

            # Mask at each height and save surface areas
            for k in np.unique(z):
                mask = z == k
                temp_profile[temp[item]][k] = area[mask].sum()
            
        results = temp_profile       
                                               
        variables[temp[item]][n, :] = results[temp[item]]

def main(t):
	vars = ('CORE_SURFACE', 'CONDENSED_SURFACE', 'PLUME_SURFACE')

	# Load cloud data for current time step
	cloud_filename = '../cloudtracker/pkl/cloud_data_%08d.pkl' % t
	clouds = cPickle.load(open(cloud_filename, 'rb'))

	ids = clouds.keys()
	ids.sort()

	data['ids'] = np.array(ids)

    # Create a savefile for surface profile at the current time step
	savefile, variables = create_savefile(t, data, vars, 'surface')

	for n, id in enumerate(ids):
		print("time: ", t, " id: ", id)
		# Select the current cloud id
		cloud = clouds[id]
		make_profiles(variables, cloud, vars, data, n)
            
        savefile.close()