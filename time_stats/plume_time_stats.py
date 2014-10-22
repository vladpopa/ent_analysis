#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset
try:
    import cPickle as pickle
except:
    import pickle
import numpy.ma
from lib.thermo import SAM
from lib.thermo import thermo
import lib.model_param as mc

def plume_time_stats():
    sample_types = ('PLUME', 'EDGE', 'SHELL', 'ENV')

    stats_dict = {}

    for time_step in range(mc.nt):
        print('Time step number:', time_step)
        cluster_dict = {}
        
        # Open plume profile file at current timestep          
        nc_files = {}
        nc_files['PLUME'] = Dataset(
            '../time_profiles/cdf/plume_profile_%08d.nc' % time_step)    	
        
        # Mask by plume area at each height
        area = nc_files['PLUME'].variables['AREA'][:]
        mask = (area > 0.)
        area[~mask] = 0.
        # Plume area
        cluster_dict['AREA'] = area[mask]
        
        # Skip time step if no plume area detected
        if mask.sum() == 0:
            nc_files['PLUME'].close()
            stats_dict[time_step] = {}
            continue

        # Masks for cloud top and bottom
        mask_top = mask.copy()
        # Account for case where plume starts at bottom of domain
        mask_top[:, 0] = False
        mask_top[:, 1:-1] = mask[:, 1:-1] & ~mask[:, 2:] & mask[:, :-2]
        mask_top
        mask_bottom = mask.copy()
        # Account for case where plume ends at top of domain
        mask_bottom[-1] = False
        mask_bottom[:, 1:-1] = mask[:, 1:-1] & mask[:, 2:] & ~mask[:, :-2]
            
        # Open edge, shell, environment, surface and stat data files
        nc_files['EDGE'] = Dataset(
            '../time_profiles/cdf/plume_edge_profile_%08d.nc' % time_step)
        nc_files['SHELL'] = Dataset(
            '../time_profiles/cdf/plume_shell_profile_%08d.nc' % time_step)
        nc_files['ENV'] = Dataset(
            '../time_profiles/cdf/plume_env_profile_%08d.nc' % time_step)
        surface_file = Dataset(
            '../time_profiles/cdf/surface_profile_%08d.nc' % time_step)
        stat_file = Dataset(mc.get_stat())
       
        # Plume z
        z = nc_files['PLUME'].variables['z'][:]
        z = np.resize(z, mask.shape)
        cluster_dict['z'] = z[mask]

        # Plume thickness
        # Use masked arrays to preserve axes; if z_min == z_max, thickness = dz
        masked_z = np.ma.masked_where(area==0., z)
        depth = np.ones_like(z)*(masked_z.max(axis=1) - 
            masked_z.min(axis=1))[:, np.newaxis] + mc.dz
        cluster_dict['depth'] = depth[mask]

        # Plume relative humidity
        relh = nc_files['PLUME'].variables['RELH'][:]
        cluster_dict['RELH'] = relh[mask]

        # Plume boundaries relative humidity
        relh = nc_files['EDGE'].variables['RELH'][:]
        cluster_dict['RELH_PLUME_EDGE'] = relh[mask]
        relh = nc_files['SHELL'].variables['RELH'][:]
        cluster_dict['RELH_PLUME_SHELL'] = relh[mask]
        relh = nc_files['ENV'].variables['RELH'][:]
        cluster_dict['RELH_PLUME_ENV'] = relh[mask]
        
        # Plume scaled z
        z = z*mask
        zmax = np.ones_like(mask)*(z.max(1))[:, np.newaxis]
        z[~mask] = 1e10
        zmin = np.ones_like(mask)*(z.min(1))[:, np.newaxis]
        cluster_dict['z_scaled'] = ((z - zmin)/(zmax - zmin))[mask]

        # Plume density
        rho = nc_files['PLUME'].variables['RHO'][:]
        cluster_dict['RHO'] = rho[mask]

        # Plume mass flux
        mf = rho*area*nc_files['PLUME'].variables['W'][:]
        cluster_dict['MF'] = mf[mask]     

        # Plume thermodynamic variables
        for var in ('W', 'QT', 'THETAV', 'THETAL', 'QN'):
            for type in sample_types:
                temp = nc_files[type].variables[var][:]   
                cluster_dict[var + '_' + type] = temp[mask]
                if var != 'W':
                    temp = stat_file.variables[var][:]
                    if var == 'QT':
                        temp = temp/1000.
                    temp2 = nc_files[type].variables[var][:] - \
                        temp[time_step, :]
                    cluster_dict[var + '_' + type + '-MEAN'] = temp2[mask]
            
            # Compute plume-shell and plume-environment differences                    
            cluster_dict[var + '_PLUME-ENV'] = cluster_dict[var + '_PLUME'] - \
                cluster_dict[var + '_ENV']
            cluster_dict[var + '_PLUME-SHELL'] = cluster_dict[var + '_PLUME'] - \
                cluster_dict[var + '_SHELL']

        # d theta_v/dz mean
        tv = stat_file.variables['THETAV'][time_step, :]
        tv[1:-1] = (tv[2:]-tv[:-2])/mc.dz/2.
        tv = tv*np.ones_like(mf) 
        cluster_dict['dTHETAV_dz_MEAN'] = tv[mask]
        
        # Plume dw/dz, dp/dz, d theta_v/dz
        for var in ('DWDZ', 'DPDZ', 'THETAV_LAPSE'):
            temp = nc_files['PLUME'].variables[var][:]
            cluster_dict[var + '_PLUME'] = temp[mask]

        # Plume surface
        surface = surface_file.variables['PLUME_SURFACE'][:]
        cluster_dict['SURFACE'] = surface[mask]

        # Height derivatives of mass flux and area
        for var in (('MF', mf), ('AREA', area)):
            temp = var[1]
            # Leapfrog, forward Euler for top and bottom 
            temp_result = (temp[:, 2:] - temp[:, :-2])/mc.dz/2.
            temp_top = (temp[:, 2:] - temp[:, 1:-1])/mc.dz
            temp_bottom = (temp[:, 1:-1] - temp[:, :-2])/mc.dz
            temp_result[mask_top] = temp_top[mask_top]
            temp_result[mask_bottom] = temp_bottom[mask_bottom]
            cluster_dict['d_' + var[0] + '_dz'] = temp_result[mask]

        # Plume time
        cluster_dict['TIME'] = np.ones_like(z[mask])*time_step*mc.dt
        
        # Data formatting
        for item in cluster_dict:
            if item in stats_dict:
                stats_dict[item].append(cluster_dict[item])
            else:
                stats_dict[item] = [cluster_dict[item]]

        # Close netCDF files
        for type in sample_types:
            nc_files[type].close()
            
    # Data formatting
    for item in stats_dict:
        stats_dict[item] = np.hstack(stats_dict[item])
    
    # Save as recarray
    n = len(stats_dict['z'])
    stats = np.zeros(1, dtype = [(key, 'f8', n) for key in stats_dict])
    for key in stats_dict:
        stats[key] = stats_dict[key]

    np.save('npy/plume_time_stats.npy', stats)

if __name__ == "__main__":
    plume_time_stats()