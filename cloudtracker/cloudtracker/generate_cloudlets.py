#!/usr/bin/env python
"""
This program generates a pkl file containing a list of dictionaries.
Each dictionary in the list represents a condensedlet.
The dictionaries have the structure:
{'core': array of ints of core points,
'condensed': array of ints of condensed points,
'plume': array of ints of plume points,
'u_condensed': ,
'v_condensed': ,
'w_condensed': ,
'u_plume': ,
'v_plume': ,
'w_plume': }
pkl files are saved in pkl/ subdirectory indexed by time
"""

import numpy
import cPickle
from utility_functions import index_to_zyx, expand_indexes

#-------------------

def expand_cloudlet(cloudlet, indexes, MC):
    """
    Given an array of indexes composing a cloudlet and a boolean mask array 
    indicating if each model index may be expanded into (True) or not (False), 
    expand the cloudlet into the permissible indices that are adjacent to the 
    cloudlet.
    
    Return an array of the indices composing the expanded cloudlet and an array 
    of the remaining indices that may be expanded into.
    """

    # Expand the cloudlet indexes into their nearest neighbours; conditional
    # sampling not taken into account.
    expanded_cloudlet = expand_indexes(cloudlet, MC)

    # Find the mask values of the expanded indexes
    mask = indexes[expanded_cloudlet]

    # Select the expanded cloudlet indexes that may be expanded into
    new_points = expanded_cloudlet[mask]

    # Remove the indices that have been added to the cloudlet
    indexes[new_points] = False

    return new_points, indexes

#---------------------

def expand_current_cloudlets(key, cloudlets, mask, MC):

    # List of all initial cloudlet points to expand from, each in an array
    cloudlet_points = []
    for cloudlet in cloudlets:
        cloudlet_points.append([cloudlet[key]])

    # Index in list of cloudlets
    cloudlet_expand_indexes = range(len(cloudlet_points))

    # Iterate while expansion is possible; start with initial cloudlets
    # Continue with list of cloudlets that expanded during the previous iteration
    while cloudlet_expand_indexes:
        next_loop_cloudlet_list = []
            
        # Go through the current list of cloudlets
        for n in cloudlet_expand_indexes:
            # Expand current cloudlet
            expanded_points, mask = expand_cloudlet(
                cloudlet_points[n][-1], mask, MC)

            # Add points expanded into to the cloudlet; if we're still 
            # expanding, try expanding this cloudlet again
            if len(expanded_points) > 0:
                cloudlet_points[n].append(expanded_points)
                next_loop_cloudlet_list.append(n)
                
        cloudlet_expand_indexes = next_loop_cloudlet_list
                
    for n, cloudlet in enumerate(cloudlet_points):
        cloudlets[n][key] = numpy.hstack(cloudlet)

    return cloudlets, mask

#---------------------

def make_new_cloudlets(key, mask, MC):
    # Create array of core/cloud/plume indices
    indexes = numpy.arange(MC['nx']*MC['ny']*MC['nz'])[mask]
    # Initialize list of cloudlets
    cloudlets = []

    # Starting at each core/cloud/plume point
    for n in indexes:
        # Clear starting point in mask
        if mask[n]:
            mask[n] = False
            
            # Start cloudlet index array with initial index
            cloudlet_indexes = [numpy.array((n,))]
            
            # Add new cloudlet; recursion until no new points in cloudlet
            done = False            
            while not done:
                # Expand cloudlet from current index
                new_indexes, mask = expand_cloudlet(cloudlet_indexes[-1], mask, MC)

                # Keep track of expanded indices in cloudlet index array
                if len(new_indexes) > 0:
                    cloudlet_indexes.append(new_indexes)
                else:
                    # Done if the number of points in the cloudlet is unchanged
                    done = True
            
            # Add a cloudlet dictionary for the new cloudlet
            cloudlet = {}
            cloudlet[key] = numpy.hstack(cloudlet_indexes)
            # Add cloudlet dictionary to the list of cloudlets
            cloudlets.append(cloudlet)
            
    return cloudlets

#-----------------

def find_mean_cloudlet_velocity(cloudlets, u, v, w, MC):
    """Calculate mean velocities of condensed region and dry plume region for 
    each cloudlet.
    """
    
    dx, dy, dz, dt = MC['dx'], MC['dy'], MC['dz'], MC['dt']                                
    ug, vg = MC['ug'], MC['vg']
    
    for cloudlet in cloudlets:
        # Condensed region
        if len(cloudlet['condensed']) > 0:
            K, J, I = index_to_zyx(cloudlet['condensed'], MC)
            # Find the mean motion of the cloudlet; account for geostrophic
            # translation of domain in SAM
            u_mean = u[K, J, I].mean()-ug
            v_mean = v[K, J, I].mean()-vg
            w_mean = w[K, J, I].mean()
            # Units are grid cells/snapshot interval
            cloudlet['u_condensed'] = round(u_mean*dt/dx)
            cloudlet['v_condensed'] = round(v_mean*dt/dy)
            cloudlet['w_condensed'] = round(w_mean*dt/dz)
        # No condensed region
        else:
            cloudlet['u_condensed'] = 0.
            cloudlet['v_condensed'] = 0.
            cloudlet['w_condensed'] = 0.

        # Plume region
        K, J, I = index_to_zyx(cloudlet['plume'], MC)
        # Find the mean motion of the cloudlet; account for geostrophic
        # translation of domain in SAM
        u_mean = u[K, J, I].mean()-ug
        v_mean = v[K, J, I].mean()-vg
        w_mean = w[K, J, I].mean()
        # Units are grid cells/snapshot interval        
        cloudlet['u_plume'] = round(u_mean*dt/dx)
        cloudlet['v_plume'] = round(v_mean*dt/dy)
        cloudlet['w_plume'] = round(w_mean*dt/dz)

    return cloudlets

#----------------------------

def generate_cloudlets(core, condensed, plume, u, v, w, MC): 
    # Find the indexes of all the core, condensed and plume points
    core = core.flatten()
    condensed = condensed.flatten()
    plume = plume.flatten()

    # Conditionally sampled fields shouldn't overlap to start with
    # Plume doesn't overlap condensed
    plume[condensed] = False
    # Condensed doesn't overlap core
    condensed[core] = False
    
    # Create the list that will hold the cloudlets; add core cloudlets
    cloudlets = make_new_cloudlets('core', core, MC)
                    
    # Initialize condensed region of each cloudlet with the core region
    for cloudlet in cloudlets:
        cloudlet['condensed'] = cloudlet['core'][:]
            
    # All current cloudlets have a core
    ncore = len(cloudlets)
    print "\t%d core cloudlets" % ncore
    
    # Expand condensed region of cloudlets that have a core
    cloudlets, condensed = expand_current_cloudlets(
        'condensed', cloudlets, condensed, MC)

    # Add any remaining points that have not been added to cloudlets 
    # as new cloudlets i.e. condensed cloudlets with no core points are new 
    # cloudlets
    condensed_cloudlets = make_new_cloudlets('condensed', condensed, MC)

    # Add condensed only cloudlets to list; these have no core indices
    for cloudlet in condensed_cloudlets:
        cloudlet['core'] = numpy.array([], dtype=numpy.int)
        cloudlets.append(cloudlet)

    # Initialize plume region of each cloudlet with the condensed region
    for cloudlet in cloudlets:
        cloudlet['plume'] = cloudlet['condensed'][:]

    # Calculate number of condensed cloudlets with no core
    ncondensed = len(cloudlets)
    print "\t%d condensed cloudlets" % (ncondensed-ncore)

    # Expand plume region of cloudlets that have a condensed region
    cloudlets, plume = expand_current_cloudlets('plume', cloudlets, plume, MC)

    # Add any remaining points that have not been added to cloudlets 
    # as new cloudlets i.e. plume cloudlets with no condensed points are new 
    # cloudlets
    plume_cloudlets = make_new_cloudlets('plume', plume, MC)

    # Add plume only cloudlets to list; these have no core or condensed indices
    for cloudlet in plume_cloudlets:
        cloudlet['core'] = numpy.array([], dtype=numpy.int)
        cloudlet['condensed'] = numpy.array([], dtype=numpy.int)
        cloudlets.append(cloudlet)

    # Calculate number of plume only cloudlets
    nplume = len(cloudlets) 
    print "\t%d plume cloudlets" % (nplume-ncondensed)

    # Calculate mean velocities of condensed region and dry plume region for 
    # each cloudlet
    cloudlets = find_mean_cloudlet_velocity(cloudlets, u, v, w, MC)

    return cloudlets

if __name__ == "__main__":
    import doctest
    doctest.testmod()