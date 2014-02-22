#!/usr/bin/env python

from __future__ import division, print_function
import sys
sys.path.append('../lib')
import numpy as np
import netCDF4
from scipy import polyfit
import model_param as mc
import glob
import matplotlib.pyplot as plt
import os
       
if __name__ == "__main__":
    snapshot_areas = np.load('npy/snapshot_shaded_area.npy')  
    snapshot_l = np.sqrt(snapshot_areas*mc.dx*mc.dy)
    
    tracked_areas = np.load('npy/tracked_shaded_area.npy')  
    tracked_l = np.sqrt(tracked_areas*mc.dx*mc.dy)
    
    bins = np.linspace(0., 4000., 100)

    # Remove extrema
    # mask1 = (tracked_l >= tracked_l.min()) & (tracked_l <= tracked_l.max())

    snapshot_h, x = np.histogram(snapshot_l, bins=bins, normed=False)
    tracked_h, x = np.histogram(tracked_l, bins=bins, normed=False)

    # tracked_h = tracked_h*1./tracked_h.sum()
    # snapshot_h = snapshot_h*1./snapshot_h.sum()

    bin_x = (x[1:] + x[:-1])/2.

    mask = (bin_x > 100) & (bin_x < 900)
    snapshot_slope, snapshot_intercept = polyfit(np.log(bin_x[mask]), \
        np.log(snapshot_h[mask]), 1)
    print('Snapshot shaded area histogram slope: {0:.2f}' \
        .format(snapshot_slope))
        
    mask = (bin_x > 200) & (bin_x < 900)
    tracked_slope, tracked_intercept = polyfit(np.log(bin_x[mask]), \
        np.log(tracked_h[mask]), 1)
    print('Tracked shaded area histogram slope: {0:.2f}' \
        .format(tracked_slope))

    xx = [x[0]]
    snapshot_hh = [snapshot_h[0]]
    tracked_hh = [snapshot_h[0]]
    for n, item in enumerate(x[1:-1]):
        xx.append(x[n+1])
        snapshot_hh.append(snapshot_h[n])
        tracked_hh.append(tracked_h[n])
        xx.append(x[n+1])
        snapshot_hh.append(snapshot_h[n+1])
        tracked_hh.append(tracked_h[n+1])

    x = np.array(xx)
    snapshot_h = np.array(snapshot_hh)
    tracked_h = np.array(tracked_hh)
    
    fig = plt.figure(1, figsize=(3.125, 3.125))
    fig.clf()
    ax = fig.add_subplot(111)
    ax.plot(x, snapshot_h, 'k-', color='k')    
    ax.plot(bin_x, (np.exp(snapshot_intercept))*bin_x**snapshot_slope, 'r-', \
        dashes=(5,2))
    ax.plot(x, tracked_h, 'k-', color='0.5')    
    ax.plot(bin_x, (np.exp(tracked_intercept))*bin_x**tracked_slope, 'b-', \
        dashes=(2,2))
        
    l = plt.legend(('snapshots', 'linear fit, slope = -2.45', \
        'tracking algorithm', 'linear fit, slope = -X.XX'), 3, prop={'size':8})
    l.draw_frame(False)
    
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$a^{1/2}$ (m)', fontsize=8)
    ax.set_ylabel('Number of Clouds per 10 m Bin', fontsize=8)
    ax.set_title('Cloud Size Distribution, SST = 297 K', fontsize=8)
    
    ax.set_xticks((50, 100, 500, 1000, 2000))
    ax.set_xticklabels(('50', '100', '500', '1000', '2000'))
    ax.set_yticks((1e1, 1e2, 1e3, 1e4))    

    ax.set_xlim(20, 4000)
    ax.set_ylim(1, 2e3)    

    # Save plot in png and eps formats
    if not os.path.exists('png'):
        os.makedirs('png')
    plt.savefig('png/SST_297_snapshot_shaded_area_histogram.png', dpi=150)
    if not os.path.exists('eps'):
        os.makedirs('eps')    
    plt.savefig('eps/SST_297_snapshot_shaded_area_histogram.eps')      
     
    plt.show()