#!/usr/bin/env python
from pylab import *
import numpy as np
import glob
from netCDF4 import Dataset
try:
    import cPickle as pickle
except:
    import pickle
import numpy.ma
from ent_analysis.lib.thermo import SAM
from ent_analysis.lib.thermo import thermo
import ent_analysis.lib.model_param as mc

def main():
    sample_types = ('CONDENSED', 'EDGE', 'SHELL', 'ENV',)

    stats_dict = {}

    for l in range(mc.nt):
        print l
        cluster_dict = {}          
        nc_files = {}
        nc_files['CONDENSED'] = Dataset('../time_profiles/cdf/condensed_profile_%08d.nc' % l)
        nc_files['PLUME'] = Dataset('../time_profiles/cdf/plume_profile_%08d.nc' % l)    	
        area = nc_files['CONDENSED'].variables['AREA'][:]
        mask = (area > 0.)
        area[~mask] = 0.
        cluster_dict['AREA'] = area[mask]
        
        if mask.sum() == 0:
            nc_files['CONDENSED'].close()
            stats_dict[l] = {}
            continue

        mask_top = mask.copy()
        mask_top[:, 1:-1] = mask[:, 1:-1] & ~mask[:, 2:] & mask[:, :-2]
        mask_bottom = mask.copy()
        mask_bottom[:, 1:-1] = mask[:, 1:-1] & mask[:, 2:] & ~mask[:, :-2]
        nc_files['EDGE'] = Dataset('../time_profiles/cdf/condensed_edge_profile_%08d.nc' % l)
        nc_files['SHELL'] = Dataset('../time_profiles/cdf/condensed_shell_profile_%08d.nc' % l)
        nc_files['ENV'] = Dataset('../time_profiles/cdf/condensed_env_profile_%08d.nc' % l)
        entrain_file = Dataset('../time_profiles/cdf/condensed_entrain_profile_%08d.nc' % l)
        surface_file = Dataset('../time_profiles/cdf/surface_profile_%08d.nc' % l)
        condensed_shell_file = Dataset('../time_profiles/cdf/condensed_shell_profile_%08d.nc' % l)
        chi_file = Dataset('../time_profiles/cdf/condensed_chi_profile_%08d.nc' % l)
        stat_file = Dataset(mc.get_stat())
       
        z = nc_files['CONDENSED'].variables['z'][:]
        z = np.resize(z, mask.shape)
        cluster_dict['z'] = z[mask]

        # Calculate and store cloud thickness for each sample
        # Use masked arrays to preserve axes; if z_min == z_max, thickness = dz
        masked_z = np.ma.masked_where(area==0., z)
        depth = np.ones_like(z)*(masked_z.max(axis=1) - 
            masked_z.min(axis=1))[:, np.newaxis] + mc.dz
        cluster_dict['depth'] = depth[mask]

        # Condensed relative humidity
        relh = nc_files['CONDENSED'].variables['RELH'][:]
        cluster_dict['RELH'] = relh[mask]

        # Condensed boundaries relative humidity
        relh = nc_files['EDGE'].variables['RELH'][:]
        cluster_dict['RELH_CONDENSED_EDGE'] = relh[mask]
        relh = nc_files['SHELL'].variables['RELH'][:]
        cluster_dict['RELH_CONDENSED_SHELL'] = relh[mask]
        relh = nc_files['ENV'].variables['RELH'][:]
        cluster_dict['RELH_CONDENSED_ENV'] = relh[mask]

        z = z*mask      
        zmax = np.ones_like(mask)*(z.max(1))[:, np.newaxis]        
        z[~mask] = 1e8
        zmin = np.ones_like(mask)*(z.min(1))[:, np.newaxis]
        cluster_dict['z_scaled'] = ((z - zmin.min())/(zmax-zmin.min()))[mask]

        rho = nc_files['CONDENSED'].variables['RHO'][:]
        cluster_dict['RHO'] = rho[mask]
        
        mf = rho*area*nc_files['CONDENSED'].variables['W'][:]
        cluster_dict['MF'] = mf[mask]

        for var in ('W', 'QT', 'THETAV', 'THETAL', 'QN'):
            for type in sample_types:
                temp = nc_files[type].variables[var][:]
        
                cluster_dict[var + '_' + type] = temp[mask]
                if var != 'W':
                    temp = stat_file.variables[var][:]
                    if var == 'QT': 
                        temp = temp/1000.
                    temp2 = nc_files['CONDENSED'].variables[var][:] - temp[l, :]
                    cluster_dict[var + '_' + type + '-MEAN'] = temp2[mask]
                                
            cluster_dict[var + '_CONDENSED-ENV'] = cluster_dict[var + '_CONDENSED'] - cluster_dict[var + '_ENV']
            cluster_dict[var + '_CONDENSED-SHELL'] = cluster_dict[var + '_CONDENSED'] - cluster_dict[var + '_SHELL']

        tv = stat_file.variables['THETAV'][l, :]
        tv[1:-1] = (tv[2:]-tv[:-2])/mc.dz/2.
        tv = tv*ones_like(temp)
        cluster_dict['dTHETAV_dz_MEAN'] = tv[mask]

        for var in ('DWDZ', 'DPDZ', 'THETAV_LAPSE'):
            temp = nc_files['CONDENSED'].variables[var][:]
            cluster_dict[var + '_CONDENSED'] = temp[mask]
            
        chi = chi_file.variables['chi'][:]
        # chi_shell = chi_file.variables['chi_shell'][:]
        # chi_se = chi_file.variables['chi_se'][:]

        cluster_dict['CHI'] = chi[mask]
        # cluster_dict['CHI_SHELL'] = chi_shell[mask]

        surface = surface_file.variables['CONDENSED_SURFACE'][:]
        cluster_dict['SURFACE'] = surface[mask]

        lsmf = stat_file.variables['MFTETCLD'][l, :]
        lsrhoa = stat_file.variables['RHO'][l, :]*stat_file.variables['VTETCLD'][l,:]

        E = entrain_file.variables['ETETCLD'][:]
        D = entrain_file.variables['DTETCLD'][:]
        massflux = entrain_file.variables['MFTETCLD'][:]
        volume = entrain_file.variables['VTETCLD'][:]
        epsilon = E/massflux
        delta = D/massflux
        wepsilon = E/rho/volume
        wdelta = D/rho/volume
               
        cluster_dict['E'] = E[mask]
        cluster_dict['D'] = D[mask]
        cluster_dict['EPSILON'] = epsilon[mask]
        cluster_dict['DELTA'] = delta[mask]
        cluster_dict['EPSILON_LS'] = (E/lsmf)[mask]
        cluster_dict['DELTA_LS'] = (D/lsmf)[mask]
        cluster_dict['WEPSILON'] = wepsilon[mask]
        cluster_dict['WDELTA'] = wdelta[mask]
        cluster_dict['WEPSILON_LS'] = (E/lsrhoa)[mask]
        cluster_dict['WDELTA_LS'] = (D/lsrhoa)[mask]

        for var in (('MF', mf), ('AREA', area)):
            temp = var[1]
            temp_result = (temp[:, 2:] - temp[:, :-2])/mc.dz/2.
            temp_top = (temp[:, 2:] - temp[:, 1:-1])/mc.dz
            temp_bottom = (temp[:, 1:-1] - temp[:, :-2])/mc.dz
            temp_result[mask_top] = temp_top[mask_top]
            temp_result[mask_bottom] = temp_bottom[mask_bottom]
            cluster_dict['d_' + var[0] + '_dz'] = temp_result[mask]

        cluster_dict['TIME'] = ones_like(z[mask])*l*mc.dt
        
        for item in cluster_dict:
            if item in stats_dict:
                stats_dict[item].append(cluster_dict[item])
            else:
                stats_dict[item] = [cluster_dict[item]]

        for type in sample_types:
            nc_files[type].close()
        entrain_file.close()
        chi_file.close()
        
    for item in stats_dict:
        stats_dict[item] = np.hstack(stats_dict[item])

    # Save as recarray
    n = len(stats_dict['z'])
    stats = np.zeros(1, dtype = [(key, 'f8', n) for key in stats_dict])
    for key in stats_dict:
       stats[key] = stats_dict[key]
    np.save('npy/condensed_time_stats.npy', stats)
        
if __name__ == "__main__":
    main()

