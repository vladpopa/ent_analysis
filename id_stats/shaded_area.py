#!/usr/bin/env python

"""
Calculate cloud shaded area for individual clouds. Following Dawe and Austin, 
2012 the functions take the projection of the clouds' condensed regions onto a
horizontal plane. Separate functions are provided for areas produced by the SAM
cloud tracking algorithm and using cloud condensed region areas taken directly
from snapshots of SAM model output.

Examples
-------
$ python shaded_area.py
$ python shaded_area.py --tracked
$ python shaded_area.py -t
"""

from __future__ import division, print_function
import sys
sys.path.append('../lib')
import glob
import os
import numpy as np
import numpy.ma as ma
import model_param as mc
import argparse
from netCDF4 import Dataset
import scipy.ndimage.measurements as measurements
import scipy.ndimage.morphology as morphology
try:
    import cPickle as pickle
except:
    import pickle
       
def snapshot_shaded_area(file_name):
    """
    Return an array of shaded cloud areas for a 3D model snapshot.
        
    Parameters
    ----------
    file_name : netCDF file name
        Data has dimensions float x(x), float y(y) and float z(z) and the field
        has been conditionally sampled for condensed area in variable int
        condensed(z, y, x).
      
    Return
    ------
    numpy.array
        Array of shaded cloud areas in units of grid cell.    
    """
    
    # Read netCDF dataset
    data = Dataset(file_name)
    condensed = data.variables['condensed'][...]
    
    # Determine shaded area mask (logical or along vertical axis)
    shaded_mask = np.any(condensed, axis=0)
    
    # Tag clouds; consider clouds connected even if they touch diagonally
    s = morphology.generate_binary_structure(2, 2)
    clouds, n = measurements.label(shaded_mask, s)
    
    # Count number of grid cells in each cloud; remove cloud free area
    cloud_areas = np.bincount(clouds.flatten().astype(int))[1:]
    
    return cloud_areas

def tracked_shaded_area(file_name):
    """
    Return an array of shaded cloud areas for a clouds tracked using the SAM
    cloud tracking algorithm.
        
    Parameters
    ----------
    file_name : pkl file name
        Data contains dictionary with cloud id as key. Values are dictionaries
        with profile type as key ('core', 'condensed', etc.) and index of grid
        points for the profile as values. Cloud id -1 represents cloud noise.
      
    Return
    ------
    numpy.array
        Array of shaded cloud areas in units of grid cell.
    """

    # Read pkl dataset 
    clouds = pickle.load(open(file_name, 'rb'))
    
    # Iterate through cloud id, calculate shaded area for each cloud id
    areas = {}
    for id, profile in clouds.iteritems():
        ij = clouds[id]['condensed']%(mc.nx*mc.ny)
        areas[id] = len(np.unique(ij))
    
    # Remove cloud noise and clouds with no condensed points
    areas.pop(-1, None)
    cloud_areas = np.ma.masked_equal(areas.values(), 0.)
      
    return cloud_areas.compressed()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine cloud shaded area')
    parser.add_argument('-t', '--tracked', help='use tracked clouds', \
        action="store_true")  
    args = parser.parse_args()
      
    # Create npy directory to store cloud areas
    if not os.path.exists('npy'):
        os.makedirs('npy')
    
    # Collect all files containing conditionally sampled condensed points
    if args.tracked:
        input_files = os.path.join('../analysis/bomex/cloudtracker/pkl', 'cloud_data_*.pkl')   
        file_list = glob.glob(input_files)
        file_list.sort()
    else:
        input_files = os.path.join(mc.data_directory, 'tracking', \
            'cloudtracker_input_*.nc')
        file_list = glob.glob(input_files)
        file_list.sort()
    
    # Calculate and store cloud shaded areas for all clouds
    cloud_areas = {}
    for n, file_name in enumerate(file_list):
        if args.tracked:
			print("Calculating shaded areas at step %d" % int(file_name[-11:-4]))
			cloud_areas[n] = tracked_shaded_area(file_name)
        else:
			print("Calculating shaded areas at step %d" % int(file_name[-11:-3]))
			cloud_areas[n] = snapshot_shaded_area(file_name)

    # Save cloud shaded areas
    if args.tracked:
        np.save('npy/tracked_shaded_area.npy', np.hstack(cloud_areas.values()))
    else:
        np.save('npy/snapshot_shaded_area.npy', np.hstack(cloud_areas.values()))