#!/usr/bin/env python

from pylab import *
import numpy as np
import scipy.stats

# import BOMEX
# import GCSSARM

def main():
    # id_stat_file = open('%s/ID_STATS/pkl/stats.pkl' % BOMEX.analysis_dir, 'rb')
    # id_stats = cPickle.load(id_stat_file)
    # cloud_stat = np.load('npy/cloud_statistics_lifetimes.npy')
    # print cloud_stat.keys()

    #     print id_stats['l_min'].min()
    # mask = (id_stats['l_min'] > 0) & (id_stats['l_max'] < (BOMEX.nt-1))
    # 
    # print mask.sum()
    # print len(mask)
    # 
    # height = id_stats['max_height'][mask]
    # base = id_stats['min_height'][mask]
    # mass = id_stats['average_mass'][mask]
    # lifetime = id_stats['duration'][mask]

    height = np.load('npy/cloud_statistics_tops.npy')
    base = np.load('npy/cloud_statistics_bases.npy')
    lifetime = np.load('npy/cloud_statistics_lifetimes.npy')
    mass = np.load('npy/cloud_statistics_masses.npy')

    print "height:", mean(height)
    print "base:", mean(base)
    print "lifetime:", mean(lifetime)/60.

    print sum(height > 1500)
    print sum(height > 1000)
    print sum(base < 1000)
    print sum(mass < 1e6)
    print sum(mass > 3e7)
    print "mass skew:", scipy.stats.skew(mass)
    print "lifetime skew:", scipy.stats.skew(lifetime)
    print "height skew:", scipy.stats.skew(height)
    print "base skew:", scipy.stats.skew(base)

    figure(figsize=(3.25, 3.5))

    subplots_adjust(left=.11, right=.99, bottom=.09, top=.95, wspace=.15, hspace=.4)

    subplot(221)
    ybins = linspace(0, 300, 31)
    
    hist(lifetime/60., bins=ybins, log=True, ec='none', color='#105ba4')
    xlabel('lifetime (min)', fontsize=6)
    ylabel('number of clouds per bin width', fontsize=6)
    title('a) Cloud Lifetime', fontsize=8)
    xticks((0, 15, 30, 45, 60, 100, 300), ('0', '15', '30', '45', '60', '100', '300'), fontsize=6)
    yticks((1, 10, 100, 1000), ('1', '10', '10$^2$', '10$^3$'), fontsize=6)
    ylim(.5, 3000)
    

    subplot(223)
    ybins = linspace(500, 3000, 31)
    
    hist(height, bins=ybins, log=True, ec='none', color='#105ba4')
    xlabel('height (km)', fontsize=6)
    ylabel('number of clouds per bin width', fontsize=6)
    title('c) Maximum Cloud Height', fontsize=8)
    xticks((500, 1000, 1500, 2000, 2500, 3000), ('.5', '1', '1.5', '2', '2.5', '3'), fontsize=6)
    yticks((1, 10, 100, 1000), ('1', '10', '10$^2$', '10$^3$'), fontsize=6)
    ylim(.5, 3000)
    xlim(500, 3000)

    subplot(224)

    ybins = linspace(500, 3000, 31)
    
    hist(base, bins=ybins, log=True, ec='none', color='#105ba4')
    xlabel('height (km)', fontsize=6)
#    ylabel('Number of Clouds', fontsize=8)
    title('d) Minimum Cloud Base', fontsize=8)
    xticks((500, 1000, 1500, 2000, 2500, 3000), ('.5', '1', '1.5', '2', '2.5', '3'), fontsize=6)
    yticks((1, 10, 100, 1000), ('', '', '', ''), fontsize=6)
    xlim(500, 3000)
    ylim(.5, 3000)

    subplot(222)

    ybins = linspace(0, 3e7, 31)
    #ybins = 50
    
    hist(mass, bins=ybins, log=True, ec='none', color='#105ba4')
    xlabel('mass (10$^7$ kg)', fontsize=6)
#    ylabel('Number of Clouds', fontsize=8)
    title('b) Mean Cloud Mass', fontsize=8)
    xticks((0, 1e7, 2e7, 3e7), ('0', '1', '2', '3'), fontsize=6)
    yticks((1, 10, 100, 1000), ('', '', '', ''), fontsize=6)
    xlim(0, 3e7)
    ylim(.5, 3000)


    plt.show()
    
    savefig('png/cloud_stats.png', dpi=300)
    savefig('eps/cloud_stats.eps')
    
if __name__ == "__main__":
    main()
