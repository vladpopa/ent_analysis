#!/usr/bin/env python

import sys, glob, os
from pylab import *
import numpy
import cPickle
try:
	from netCDF4 import Dataset
except:
	try:
		from netCDF3 import Dataset
	except:
		from pupynere import netcdf_file as Dataset
import networkx

import ent_analysis.lib.model_param as mc

def main(item):
    created_file_ids = []
    for t in range(mc.nt):
        ncfile = Dataset('../time_profiles/cdf/%s_profile_%08d.nc' % (item, t))
        print 'time_profiles/cdf/%s_profile_%08d.nc' % (item, t)

        ids = ncfile.variables['ids'][:]
        z = ncfile.variables['z'][:]

        for n, id in enumerate(ids):
            id = int(id)
            print "time: ", t, " id: ", id
            if id not in created_file_ids:
                savefile = Dataset('cdf/%s_profile_%08d.nc' % (item, id), 'w')
                
                # Create savefile
                savefile.createDimension('t', None)
                savefile.createDimension('z', len(z))

                tsavevar = savefile.createVariable('t', 'd', ('t',))
                tsavevar[0] = t
                zsavevar = savefile.createVariable('z', 'd', ('z',))
                zsavevar[:] = z[:]

                for name in ncfile.variables:
                    if name not in ('ids', 'z'):
                        new_variable = savefile.createVariable(name, 'd', ('t', 'z'))
                        new_variable[0, :] = ncfile.variables[name][n, :]
                savefile.close()
                created_file_ids.append(id)
            else:
                savefile = Dataset('cdf/%s_profile_%08d.nc' % (item, id), 'a')
                tvar = savefile.variables['t']
                l = len(tvar)
                tvar[l] = t
                for name in ncfile.variables:
                    if name not in ('ids', 'z'):
                        savefile.variables[name][l, :] = ncfile.variables[name][n, :]
                savefile.close()
                  
        ncfile.close()
   
if __name__ == "__main__":
    for item in ('condensed', 'condensed_env', 
      'condensed_edge', 'condensed_shell',
      'core','core_env', 
      'core_edge', 'core_shell',
      'plume', 'surface', 'condensed_entrain', 'core_entrain'):
        main(item)
