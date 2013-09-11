import numpy as np
import site
import sys
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
import cPickle
import glob
from cloudtracker.utility_functions import find_halo
from util_functions import index_to_xy, label_cloud_depths, find_cloud_indices, label_depth, \
                           find_cloud_ids_at_z, find_indices_at_z, \
                           label_cloud_depths_at_z, filter_clusters


"""

returns cloud lifetimes, areas

inputs: filenames  - a list of .pkl files, each file contains a dictionary of Cluster objects (from cloudtracker output)
                     note: file names should be time stamped
        filtered_ids - ids of clusters after filtering
        maxid - maximum cluster id
        MC - a dictionary with the model configuration
    
                   
"""

def compute_cloud_vars(filenames, filtered_ids, maxid, MC):

    filenames.sort()

    #lifetimes of clouds
    cloud_lifetimes = np.zeros(maxid+1)
    #areas of xy projections of clouds
    areas = np.zeros(maxid+1)

    for fname in filenames:
        clusters = cPickle.load(open(fname, 'rb'))
        #iterate over each cluster at the current time
        #increment the lifetime if it has a condensed region,
        #and find the area of xy cloud projection
        for id, cluster in clusters.iteritems():
            if id in filtered_ids and cluster.has_condensed():
                cloud_lifetimes[id] = cloud_lifetimes[id] + MC['dt']
                proj = index_to_xy(cluster.condensed_mask(), MC).T
                areas[id] = areas[id] + len(proj)

    cloud_ids = np.where(cloud_lifetimes > 0)[0]
    #cloud_ids = np.intersect1d(cloud_ids, filtered_ids)

    cloud_lifetimes = cloud_lifetimes[cloud_ids]
    areas = areas[cloud_ids]
    #average the cloud projection areas over their lifetime
    areas = areas*(MC['dx']*MC['dy'])/(cloud_lifetimes/MC['dt'])

    return cloud_lifetimes, areas


"""

 inputs: filename  -  name of .pkl file at some time step (contains a dictionary of Cluster objects)
         filtered_ids - ids of clusters after filtering
         MC - a dictionary with the model configuration
         h (optional parameter) - height level
        

 output: an array of distances to cloud edge (the index gives the grid point index)

"""

def compute_distances_to_cloud_edges(filename, filtered_ids, MC, h=None):
   

    clusters = cPickle.load(open(filename, 'rb'))
    
    cloud_depths = label_cloud_depths(clusters, filtered_ids, MC)
    cloud_indices, env_indices = find_cloud_indices(clusters, filtered_ids, MC)
    env_halo = find_halo(env_indices, MC)
    env_depth = label_depth(env_indices, env_halo, MC)

    if h == None:
        #label cloud/environment distances to cloud edge
        distances = np.zeros(MC['ny']*MC['nx']*MC['nz'])
        for id, depths in cloud_depths.iteritems():
            for n, indices in enumerate(depths):
            #assuming grid spacing is equal in all directions
                distances[indices] = -(n+.5)*MC['dx']
            for n, indices in enumerate(env_depth):
                distances[indices] = (n+.5)*MC['dx']
            return distances
    else:
    #label cloud/environment distances to cloud edge at height level h
        cloud_ids_at_h = find_cloud_ids_at_z(h, clusters, MC)
        cloud_ids_at_h = np.intersect1d(filtered_ids, cloud_ids_at_h)
        #if there are clouds at height h, compute distances to cloud edges
        #if len(cloud_ids_at_h):
        cloud_depths_at_h = label_cloud_depths_at_z(h, clusters, cloud_ids_at_h, MC)
        env_indices_at_h = find_indices_at_z(h, env_indices, MC)
        env_halo_at_h = find_indices_at_z(h, env_halo, MC)
        env_depth_at_h = label_depth(env_indices_at_h, env_halo_at_h, MC)

        distances_at_h = np.empty(MC['ny']*MC['nx']*MC['nz'])
        distances_at_h[:] = np.NAN


        if len(cloud_ids_at_h):
            for id, depths in cloud_depths_at_h.iteritems():
                for n, indices in enumerate(depths):
                    distances_at_h[indices] = -(n+.5)*MC['dx']

                for n, indices in enumerate(env_depth_at_h):
                    distances_at_h[indices] = (n+.5)*MC['dx']

        valid_indices = np.where(np.isfinite(distances_at_h))[0]
        distances_at_h = distances_at_h[valid_indices]

        return distances_at_h, cloud_ids_at_h














    
    
    







