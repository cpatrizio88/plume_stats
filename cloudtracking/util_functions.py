import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
import cPickle
import glob
from cloudtracker.utility_functions import index_to_zyx, zyx_to_index, find_halo, expand_indexes


#library of functions to calculate plume/cloud statistics

"""
 function to project indices onto the x-y plane

 inputs: indices - indices to project
         ny, nx - grid points in y, and x directions respectively
                      (i.e. the horizontal dimensions of the array)

 output: returns an array with the x-y projections of indices
         (first row: x coordinate, second row: y coordinate)

"""

def index_to_xy(index, MC):
    #ny = MC['ny']
    #nx = MC['nx']
    #index = index % (ny*nx)
    #y = index / nx
    #x = index % nx
    #tmp = zip(x, y)
    #return np.unique(tmp).T
    points = index_to_zyx(index, MC)
    z = points[0,:]
    y = points[1,:]
    x = points[2,:]
    tmp = zip(x, y)
    return np.unique(tmp).T

def xy_to_index(x, y, MC):
   nx = MC['nx']
   return nx*y + x

def index_to_xz(index, MC):
    points = index_to_zyx(index, MC)
    z = points[0,:]
    y = points[1,:]
    x = points[2,:]
    tmp = zip(x, z)
    return np.unique(tmp).T
    

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

inputs:  clusters - a dictionary of Cluster objects at a single time step (indexed by cloud id)
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
  returns grid point indices for all clouds with ids in cloud_ids
  and the compliment of these grid point indices (the environment)
  
"""

def find_cloud_indices(clusters, cloud_ids, MC):
    
    cloud_indices = np.array([])

    grid = np.arange(MC['nx']*MC['ny']*MC['nz'])

    for id, cluster in clusters.iteritems():
        if id in cloud_ids:
            tmp = cluster.condensed_mask()
            cloud_indices = np.concatenate([cloud_indices, tmp])

    env_indices = np.setdiff1d(grid, cloud_indices, assume_unique=True)

    return cloud_indices, env_indices

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
    points = index_to_zyx(indices, MC)
    zcoords = points[0,:]
    ycoords = points[1,:]
    xcoords = points[2,:]
    hit = np.where(zcoords == h)
    xcoords = xcoords[hit]
    ycoords = ycoords[hit]
    zcoords = zcoords[hit]

    return zyx_to_index(zcoords, ycoords, xcoords, MC)


    #nx = MC['nx']
    #ny = MC['ny']
    
    #hit = np.logical_and(h*nx*ny <= indices, indices < (h+1)*nx*ny)
    #return indices[hit]

def find_indices_at_y(y, indices, MC):
    points = index_to_zyx(indices, MC)
    zcoords = points[0,:]
    ycoords = points[1,:]
    xcoords = points[2,:]
    hit = np.where(ycoords == y)
    xcoords = xcoords[hit]
    ycoords = ycoords[hit]
    zcoords = zcoords[hit]

    return zyx_to_index(zcoords, ycoords, xcoords, MC)

"""
returns ids for clouds at height level h

inputs: h - height level
        clusters - dictionary of Cluster objects at some time 
        MC - dictionary with model configuration

"""
def find_cloud_ids_at_height(h, clusters, MC):
    ids = []
    for id, cluster in clusters.iteritems():
        indices_at_h = find_indices_at_height(h, cluster.condensed_mask(), MC)
        if len(indices_at_h):
            ids.append(id)
    return ids

"""
returns ids for any plumes (that are not entirely condensed) at North-South grid point y
"""
def find_plume_ids_at_y(y, clusters, MC):
    ids = []
    for id, cluster in clusters.iteritems():
        plume = cluster.plume_mask()
        condensed = cluster.condensed_mask()
        dry_plume = np.setdiff1d(plume, condensed)
        indices_at_y = find_indices_at_y(y, dry_plume, MC)
        if len(indices_at_y):
            ids.append(id)
    return ids

"""
returns ids for plumes at timestep t

inputs: filenames - list of .pkl files, each containing a dictionary of Cluster objects
                    (from cloudtracker output)
        t - timestep
        filtered_ids - list of ids after filtering

output: list of plume ids (any plume ids that are not all condensed points)
"""
def find_plume_ids_at_t(t, filenames, filtered_ids):
    filenames.sort()
    clusters = cPickle.load(open(filenames[t], 'rb'))
    plume_ids = []
    for id, cluster in clusters.iteritems():
        if id in filtered_ids:
           dry_plume = np.setdiff1d(cluster.plume_mask(), cluster.condensed_mask())
           if len(dry_plume):
              plume_ids.append(id)
    return plume_ids
    

"""
returns a dictionary indexed by cloud id with depth arrays at height level h

inputs: h - height level 
        clusters - dictionary of Cluster objects indexed by id
        ids - ids of clouds to label depth 
        MC - dictionary with model configuration

"""
def label_cloud_depths_at_height(h, clusters, ids, MC):
    cloud_depths = {}
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


"""

returns ids of clusters that have lifetimes greater than one
model output time step, and are greater than some size threshold


inputs: filenames - list of filenames of .pkl files from cloudtracker output
                   (each .pkl file contains a dictionary of Cluster objects at some time step)

outputs: filtered_ids - a list of cluster ids that pass above criteria
         maxid - the maximum cluster id
         MC - dictionary of model configuration
                   
"""

def filter_clusters(filenames):
    
    filenames.sort()
    cluster = cPickle.load(open(filenames[0], 'rb'))[0]
    MC = cluster.MC
    #get model parameters
    dt = MC['dt']
    nx = MC['nx'] #number of grid points in x, y, and z direction
    ny = MC['ny']
    nz = MC['nz']
    nt = MC['nt'] #number of time steps
    dx = MC['dx'] #grid spacing (m)
    dy = MC['dy']
    dz = MC['dz']
    #a list of dictionaries; each dictionary contains the set of all
    #cluster objects at a given time
    clusters_list = []
    maxid = 0

    for fname in filenames:
        clusters = cPickle.load(open(fname, 'rb'))
        clusters_list.append(clusters)
        ids = clusters.keys()
        maxid_tmp = max(ids)
        if maxid_tmp > maxid:
            maxid = maxid_tmp

    #lifetimes of clusters
    lifetimes = np.zeros(maxid+1)
    #volumes of clusters (# of grid cells that the plume occupies)
    volumes = np.zeros(maxid+1)
    volume_threshold = 1

    #clusters is a dictionary of Cluster objects
    #(see cloud_objects.py in cloudtracker)
    for clusters in clusters_list:
        #iterate over each cluster at the current time
        #increment the lifetime and volume of the cluster 
        #and find the area of xy cloud projection
        for id, cluster in clusters.iteritems():
            lifetimes[id] = lifetimes[id] + dt
            volumes[id] = len(cluster.plume_mask())

    #filter out noise (i.e. clusters that exist only for a single time step)
    filtered_ids = np.where(np.logical_and(lifetimes > dt, volumes > volume_threshold))[0]
    
    return filtered_ids, maxid, MC









