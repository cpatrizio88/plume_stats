import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
from cloudtracker.utility_functions import index_to_zyx, zyx_to_index, find_halo, expand_indexes


#library of functions to calculate plume/cloud statistics

"""
 function to project indices onto the x-y plane
 note: numpy arrays are in row-major order by default
 
 inputs: indices - indices to project
         ny, nx - grid points in y, and x directions respectively
                      (i.e. the horizontal dimensions of the array)

 output: returns an array with the x-y projections of indices
         (first row: x coordinate, second row: y coordinate)

"""

def index_to_xy(index, MC):
    ny = MC['ny']
    nx = MC['nx']
    index = index % (ny*nx)
    y = index / nx
    x = index % nx
    tmp = zip(x, y)
    return np.unique(tmp).T

def xy_to_index(x, y, MC):
   nx = MC['nx']
   return nx*y + x

"""

returns a list of arrays, each array contains the grid indices for a depth level
in a cloud/cloud environment 

based off of label_distance bomex_histogram_data.py  in /users/jdawe/figures/knut

inputs: indices - grid indices specifying the cloud/environment
        MC - dictionary containing grid parameters (specified in model_config.cfg)
        halo - grid indices adjacent to the cloud boundary

output: a list of numpy arrays with grid indices at each depth level (grid indices
        that are equidistant from the cloud boundary) in a cloud/cloud environment
        (starts at grid indices closest to the cloud boundary) 
        
"""

def label_depth(indices, halo, MC):

   mask_indexes = indices[:]

   done = False
   depth = []
   while not done:
        expanded_indexes = expand_indexes(halo, MC)

        # From the expanded indexes list, select the points in the indexes list
        # expand_index_list returns unique values so we don't have to filter duplicates
        new_indexes = np.intersect1d(expanded_indexes, mask_indexes, assume_unique=True)

        if len(new_indexes) == 0:
            done = True
        else:
            depth.append(new_indexes)

            # From the mask_indexes, remove all the new points
            mask_indexes = np.setdiff1d(mask_indexes, new_indexes, assume_unique=True)

            halo = new_indexes

   # return a list composed of arrays of the points at each depth
   return depth


"""
 returns a dictionary (indexed by cloud id) containing cloud depth arrays
 (as described in label_depth) 

inputs: clusters - a dictionary of Cluster objects at a single time step (indexed by cloud id)
         cloud_ids - a list of ids for which cloud depth arrays will be computed
                     must be a subset of the id's found in clusters
         MC - dictionary containing grid parameters (as in model_config.cfg)

 outputs: see above

"""

def label_cloud_depths(clusters, ids, MC):
    cloud_depths = {}
    #env_indices, cloud_indices = find_all_indices(clusters, MC)
    for id, cluster in clusters.iteritems():
        if id in ids:
              indices = cluster.condensed_mask()
              cloud_halo = cluster.condensed_halo()
              cloud_depths[id] = label_depth(indices, cloud_halo, MC)

    return cloud_depths
              

"""
  returns all grid indices that are not condensed points
  and all grid indices that are condensed points as numpy arrays

"""

def find_all_indices(clusters, MC):
    
    cloud_indices = np.array([])

    grid = np.arange(MC['nx']*MC['ny']*MC['nz'])

    for id, cluster in clusters.iteritems():
        tmp = cluster.condensed_mask()
        cloud_indices = np.concatenate([cloud_indices, tmp])

    env_indices = np.setdiff1d(grid, cloud_indices, assume_unique=True)

    return env_indices, cloud_indices

#calculates the distance between index1 and index2
#(note: index2 or index1 may be an array of indices, but not both)
def calc_distance(index1, index2, MC):

    point1 = index_to_zyx(index1, MC)
    
    point2 = index_to_zyx(index2, MC)
    
    # Calculate distances
    ny, nx = MC['ny'], MC['nx']
    dy, dx, dz = MC['dy'], MC['dx'], MC['dz']
        
    delta_x = (point2[2, ...] - point1[2, ...])*dx

    delta_y = (point2[1, ...] - point1[1, ...])*dy
   
    delta_z =(point2[0, ...] - point1[0, ...])*dz

    return np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)

"""
returns the grid point indices at height level h

"""

def find_indices_at_height(h, indices, MC):

    nx = MC['nx']
    ny = MC['ny']
    
    hit = np.logical_and(h*nx*ny <= indices, indices < (h+1)*nx*ny)
    return indices[hit]
"""
returns ids for clouds at height level h

"""
def find_cloud_ids_at_height(h, clusters, MC):
    ids = []
    for id, cluster in clusters.iteritems():
        indices_at_h = find_indices_at_height(h, cluster.condensed_mask(), MC)
        if len(indices_at_h):
            ids.append(id)
    return ids
    

"""
returns a dictionary indexed by cloud id with depth arrays at height level h

inputs: h - height level 
        clusters - dictionary of Cluster objects indexed by id
        MC - dictionary with model configuration

"""
def label_cloud_depths_at_height(h, clusters, MC):
    cloud_depths = {}
    ids = find_cloud_ids_at_height(h, clusters, MC)
    for id, cluster in clusters.iteritems():
        if id in ids:
              indices_at_h = find_indices_at_height(h, cluster.condensed_mask(), MC)
              cloud_halo = find_halo(indices_at_h, MC)
              cloud_halo_at_h = find_indices_at_height(h, cloud_halo, MC)
              depths = label_depth(indices_at_h, cloud_halo_at_h, MC)
              cloud_depths[id] = depths
              
    return cloud_depths

def expand_indexes_horizontal(indexes, MC):
    # Expand a given set of x-y grid indexes to include the nearest
    # neighbour points in the horizontal direction
    # indexes is an array of x-y grid indexes
    
    ny, nx =  MC['ny'], MC['nx']
                    
    I_J = index_to_xy( indexes, MC )

    stack_list = [I_J, ]
    for item in ((-1, 0), (1, 0),
                 (0, -1), (0, 1)): 
          stack_list.append(I_J + np.array(item)[:, np.newaxis] )
    
    expanded_index = np.hstack(stack_list)

    # re-entrant domain
    expanded_index[0, :] = expanded_index[0, :]%nx
    expanded_index[1, :] = expanded_index[1, :]%ny

    # convert back to indexes
    expanded_index = xy_to_index(expanded_index[0, :],
                                  expanded_index[1, :],
                                 MC)
                                  
    expanded_index = np.unique(expanded_index)
    
    return expanded_index








