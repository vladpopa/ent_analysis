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
from lib.thermo import SAM
from lib.thermo import thermo
import lib.model_param as mc

def core_time_stats():
    sample_types = ('CORE', 'ENV', 'PLUME')

    stats_dict = {}

    for l in range(mc.nt):
        print l
        cluster_dict = {}
        
        nc_files = {}
        nc_files['CORE'] = Dataset('../time_profiles/cdf/core_profile_%08d.nc' % l)
        nc_files['PLUME'] = Dataset('../time_profiles/cdf/plume_profile_%08d.nc' % l)	
        area = nc_files['CORE'].variables['AREA'][:]
        mask = (area > 0.)
        area[~mask] = 0.
        cluster_dict['AREA'] = area[mask]
        
        if mask.sum() == 0:
            nc_files['CORE'].close()
            nc_files['PLUME'].close()
            continue

        nc_files['EDGE'] = Dataset('../time_profiles/cdf/core_edge_profile_%08d.nc' % l)
        nc_files['SHELL'] = Dataset('../time_profiles/cdf/core_shell_profile_%08d.nc' % l)
        nc_files['ENV'] = Dataset('../time_profiles/cdf/core_env_profile_%08d.nc' % l)
        entrain_file = Dataset('../time_profiles/cdf/core_entrain_profile_%08d.nc' % l)
        surface_file = Dataset('../time_profiles/cdf/surface_profile_%08d.nc' % l)
        condensed_shell_file = Dataset('../time_profiles/cdf/condensed_shell_profile_%08d.nc' % l)
        chi_file = Dataset('../time_profiles/cdf/core_chi_profile_%08d.nc' % l)
        stat_file = Dataset(mc.get_stat())

        z = nc_files['CORE'].variables['z'][:]
        z = np.resize(z, mask.shape)
        cluster_dict['z'] = z[mask]
        
        # Calculate and store core thickness for each sample
        # Use masked arrays to preserve axes; if z_min == z_max, thickness = dz
        masked_z = np.ma.masked_where(area==0., z)
        depth = np.ones_like(z)*(masked_z.max(axis=1) - 
            masked_z.min(axis=1))[:, np.newaxis] + mc.dz
        cluster_dict['depth'] = depth[mask]
        
        # Core relative humidity
        relh = nc_files['CORE'].variables['RELH'][:]
        cluster_dict['RELH'] = relh[mask]

        # Core boundaries relative humidity
        relh = nc_files['EDGE'].variables['RELH'][:]
        cluster_dict['RELH_CORE_EDGE'] = relh[mask]
        relh = nc_files['SHELL'].variables['RELH'][:]
        cluster_dict['RELH_CORE_SHELL'] = relh[mask]
        relh = nc_files['ENV'].variables['RELH'][:]
        cluster_dict['RELH_CORE_ENV'] = relh[mask]

        stat_core = stat_file.variables['COR'][l, :]
        if (stat_core > 0.).any():
            k_cb = np.nonzero(stat_core > 0.)[0].min()
        else:
            k_cb = np.nonzero(area)[1].min()
            
        z_cb = ones_like(z)*z[0, k_cb]
        cluster_dict['z_cb'] = z_cb[mask]
        
        z = z*mask
        zmax = np.ones_like(mask)*(z.max(1))[:, np.newaxis]        
        z[~mask] = 1e10
        zmin = np.ones_like(mask)*(z.min(1))[:, np.newaxis]
        cluster_dict['z_cb'] = ((z - zmin)/(zmax - zmin))[mask]

        rho = nc_files['CORE'].variables['RHO'][:]
        cluster_dict['RHO'] = rho[mask]
        
        mf = rho*area*nc_files['CORE'].variables['W'][:]
        cluster_dict['MF'] = mf[mask]

        for var in ('W', 'QT', 'THETAV', 'THETAL', 'QN'):
            for type in sample_types:
                temp = nc_files[type].variables[var][:]
        
                cluster_dict[var + '_' + type] = temp[mask]
                if var != 'W':
                    temp = stat_file.variables[var][:]
                    if var == 'QT': 
                        temp = temp/1000.
                    temp2 = nc_files[type].variables[var][:] - temp[l, :]
                    cluster_dict[var + '_' + type + '-MEAN'] = temp2[mask]
                                
            cluster_dict[var + '_CORE-ENV'] = cluster_dict[var + '_CORE'] - cluster_dict[var + '_ENV']
            #cluster_dict[var + '_CORE-SHELL'] = cluster_dict[var + '_CORE'] - cluster_dict[var + '_SHELL']

        qsat = SAM.qsatw(nc_files['CORE'].variables['TABS'][:], 
                         nc_files['CORE'].variables['PRES'][:])
        cluster_dict['QSAT_CORE'] = qsat[mask]

        qsat_cb = qsat[:, k_cb]
        qsat_cb = ones_like(qsat)*qsat_cb[:, np.newaxis]
        cluster_dict['QSAT_CB'] = qsat_cb[mask]
        
        tv = stat_file.variables['THETAV'][l, :]
        tv[1:-1] = (tv[2:]-tv[:-2])/mc.dz/2.
        tv = tv*ones_like(mf)
        cluster_dict['dTHETAV_dz_MEAN'] = tv[mask]

        chi = chi_file.variables['chi_theta'][:]
        cluster_dict['CHI'] = chi[mask]

        surface = surface_file.variables['CORE_SURFACE'][:]
        cluster_dict['SURFACE'] = surface[mask]

        lsmf = stat_file.variables['MFTETCOR'][l, :]
        lsrhoa = stat_file.variables['RHO'][l, :]*stat_file.variables['VTETCOR'][l,:]
        
        qc = stat_file.variables['QTCOR'][l, :]/1000.
        qe = stat_file.variables['QTCEN'][l, :]/1000.
        tc = stat_file.variables['TLCOR'][l, :]
        te = stat_file.variables['TLCEN'][l, :]
        wc = stat_file.variables['WCOR'][l, :]
        we = stat_file.variables['WCEN'][l, :]
        
        dwdt = entrain_file.variables['DWDT'][:]
        E = entrain_file.variables['ETETCOR'][:]
        D = entrain_file.variables['DTETCOR'][:]
        Eq = entrain_file.variables['EQTETCOR'][:]
        Dq = entrain_file.variables['DQTETCOR'][:]
        Et = entrain_file.variables['ETTETCOR'][:]
        Dt = entrain_file.variables['DTTETCOR'][:]
        Ew = entrain_file.variables['EWTETCOR'][:]
        Dw = entrain_file.variables['DWTETCOR'][:]
        
        massflux = entrain_file.variables['MFTETCOR'][:]
        volume = entrain_file.variables['VTETCOR'][:]
               
        cluster_dict['DWDT'] = dwdt[mask]
        cluster_dict['E'] = E[mask]
        cluster_dict['D'] = D[mask]
        cluster_dict['EQ'] = Eq[mask]
        cluster_dict['DQ'] = Dq[mask]
        cluster_dict['ET'] = Et[mask]
        cluster_dict['DT'] = Dt[mask]
        cluster_dict['EW'] = Ew[mask]
        cluster_dict['DW'] = Dw[mask]
        
        Aq = ((Eq/E) - qe)/(qc - qe)
        Bq = (qc - (Dq/D))/(qc - qe)
        At = ((Et/E) - te)/(tc - te)
        Bt = (tc - (Dt/D))/(tc - te)
        Aw = ((Ew/E) - we)/(wc - we)
        Bw = (wc -(Dw/D))/(wc - we)
        
        cluster_dict['AQ'] = Aq[mask]
        cluster_dict['BQ'] = Bq[mask]
        cluster_dict['AT'] = At[mask]
        cluster_dict['BT'] = Bt[mask]
        cluster_dict['AW'] = Aw[mask]
        cluster_dict['BW'] = Bw[mask]
        
        Eq_T = ((qc*(E-D) - (Eq-Dq))/(qc-qe))
        Dq_T = ((qe*(E-D) - (Eq-Dq))/(qc-qe))
        Et_T = ((tc*(E-D) - (Et-Dt))/(tc-te))
        Dt_T = ((te*(E-D) - (Et-Dt))/(tc-te))
        Ew_T = ((wc*(E-D) - (Ew-Dw))/(wc-we))
        Dw_T = ((we*(E-D) - (Ew-Dw))/(wc-we))
        
        cluster_dict['EQ_T'] = Eq_T[mask]
        cluster_dict['DQ_T'] = Dq_T[mask]
        cluster_dict['ET_T'] = Et_T[mask]
        cluster_dict['DT_T'] = Dt_T[mask]
        cluster_dict['EW_T'] = Ew_T[mask]
        cluster_dict['DW_T'] = Dw_T[mask]
        
        cluster_dict['EPSILON'] = (E/massflux)[mask]
        cluster_dict['Q_EPSILON'] = (Eq/massflux)[mask]
        cluster_dict['T_EPSILON'] = (Et/massflux)[mask]
        cluster_dict['W_EPSILON'] = (Ew/massflux)[mask]
        cluster_dict['Q_EPSILON_T'] = (Eq_T/massflux)[mask]
        cluster_dict['T_EPSILON_T'] = (Et_T/massflux)[mask]
        cluster_dict['W_EPSILON_T'] = (Ew_T/massflux)[mask]
        
        cluster_dict['DELTA'] = (D/massflux)[mask]
        cluster_dict['Q_DELTA'] = (Dq/massflux)[mask]
        cluster_dict['T_DELTA'] = (Dt/massflux)[mask]
        cluster_dict['W_DELTA'] = (Dw/massflux)[mask]
        cluster_dict['Q_DELTA_T'] = (Dq_T/massflux)[mask]
        cluster_dict['T_DELTA_T'] = (Dt_T/massflux)[mask]
        cluster_dict['W_DELTA_T'] = (Dw_T/massflux)[mask]

        for var in ('DWDZ', 'DPDZ', 'THETAV_LAPSE'):
            temp = nc_files['CORE'].variables[var][:]
            cluster_dict[var + '_CORE'] = temp[mask]            

        ww_reyn = nc_files['CORE'].variables['WWREYN'][:]
        wq_reyn = nc_files['CORE'].variables['WQREYN'][:]
        
        ww_reyn = ww_reyn*rho*area
        wq_reyn = wq_reyn*rho*area
        
        ww_reyn[~mask] = 0.
        wq_reyn[~mask] = 0.

        qt_core = nc_files['CORE'].variables['QT'][:]
        qt_core[~mask] = 0.

        mask_top = mask.copy()
        mask_top[:, 1:-1] = mask[:, 1:-1] & ~mask[:, 2:] & mask[:, :-2]
        mask_bottom = mask.copy()
        # Account for case where core ends at top of domain
        mask_bottom[-1] = False
        mask_bottom[:, 1:-1] = mask[:, 1:-1] & mask[:, 2:] & ~mask[:, :-2]

        for var in (('MF', mf), 
                    ('AREA', area), 
                    ('WW', ww_reyn), 
                    ('WQ', wq_reyn),
                    ('Q', qt_core),):
            temp = var[1]
            temp_result = (temp[:, 2:] - temp[:, :-2])/mc.dz/2.
            temp_top = (temp[:, 2:] - temp[:, 1:-1])/mc.dz
            temp_bottom = (temp[:, 1:-1] - temp[:, :-2])/mc.dz
            temp_result[mask_top] = temp_top[mask_top]
            temp_result[mask_bottom] = temp_bottom[mask_bottom]
            cluster_dict['D' + var[0] + 'DZ_CORE'] = temp_result[mask]

        cluster_dict['TIME'] = ones_like(z[mask])*l*mc.dt
        ids = nc_files['CORE'].variables['ids'][:]
        cluster_dict['ID'] = (ones_like(z)*ids[:, np.newaxis])[mask]
        
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
    np.save('npy/core_stats.npy', stats)
        
if __name__ == "__main__":
    core_stats()