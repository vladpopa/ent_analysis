#!/usr/bin/env python
import numpy as np
try:
	from netCDF4 import Dataset
except:
	try:
		from netCDF3 import Dataset
	except:
		from pupynere import netcdf_file as Dataset
import lib.model_param as mc

def main(item):
    # Keep track of id encountered at previous time steps
    created_file_ids = []
    for t in range(mc.nt):
        ncfile = Dataset('../time_profiles/cdf/%s_profile_%08d.nc' % (item, t))
        print 'time_profiles/cdf/%s_profile_%08d.nc' % (item, t)

        ids = ncfile.variables['ids'][:]
        z = ncfile.variables['z'][:]

        # Iterate through ids at current time step
        for n, id in enumerate(ids):
            id = int(id)
            print "time: ", t, " id: ", id
            # New id
            if id not in created_file_ids:
                # Create new id file
                savefile = Dataset('cdf/%s_profile_%08d.nc' % (item, id), 'w')
                
                savefile.createDimension('t', None)
                savefile.createDimension('z', len(z))

                tsavevar = savefile.createVariable('t', 'd', ('t',))
                tsavevar[0] = t
                zsavevar = savefile.createVariable('z', 'd', ('z',))
                zsavevar[:] = z[:]

                for name in ncfile.variables:
                    if name not in ('ids', 'z'):
                        new_variable = savefile.createVariable(
                            name, 'd', ('t', 'z'))
                        new_variable[0, :] = ncfile.variables[name][n, :]
                savefile.close()
                created_file_ids.append(id)
            # Existing id
            else:
                # Existing id file
                savefile = Dataset('cdf/%s_profile_%08d.nc' % (item, id), 'a')
                tvar = savefile.variables['t']
                l = len(tvar)
                tvar[l] = t
                for name in ncfile.variables:
                    if name not in ('ids', 'z'):
                        savefile.variables[name][l, :] = \
                            ncfile.variables[name][n, :]
                savefile.close()
                  
        ncfile.close()
   
if __name__ == "__main__":
    for item in ('plume', 'plume_env', 'plume_edge', 'plume_shell',
        'condensed', 'condensed_env', 'condensed_edge', 'condensed_shell',
        'core','core_env', 'core_edge', 'core_shell',
        'surface', 'condensed_entrain', 'core_entrain'):
        main(item)