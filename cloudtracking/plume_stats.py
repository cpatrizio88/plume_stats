from __future__ import division
import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
import cPickle
import glob
from netCDF4 import Dataset
from util_functions import index_to_xy, filter_clusters, label_depths, \
                           find_plume_indices, label_depth, find_plume_ids_at_z, find_indices_at_z, \
                           label_depths_at_z
from cloudtracker.utility_functions import find_halo



"""
  outputs statistics for plumes stored as Cluster objects in .pkl files (cloudtracker output)

  input: filenames - contains list names of .pkl files 

  outputs: array of mean lifetime (min)
           arary of mean size (m)
           array of booleans indicating whether a plume forms a condensed region or not

           note: all plumes that exist for only a single model output time step are ignored 
"""

def compute_plume_vars(filenames, filtered_ids, maxid, MC):

    filenames.sort()

    #lifetimes of clusters
    lifetimes = np.zeros(maxid+1)
   
    #areas of xy projections of plumes (before condensed region forms)
    areas = np.zeros(maxid+1)
    #keep track of plumes that form condensed regions
    #initialize all plumes to False
    #has_condensed = np.zeros(maxid+1) < -1
    #track when plumes condense (i.e. the model output timestep when condensation occurs)
    #a NaN entry signifies that the plume has not condensed
    condensed_time = np.empty(maxid+1)
    condensed_time[:] = np.NAN
    #keep track of height of condensation
    condensed_height = np.empty(maxid+1)
    condensed_height[:] = np.NAN
  
    for t, fname in enumerate(filenames):
        clusters = cPickle.load(open(fname, 'rb'))
        #iterate over each cluster at the current time
        #increment the lifetime of the plume
        #and find the area of plume x-y projection
        for id, cluster in clusters.iteritems():
           if id in filtered_ids:
               #find height and time (timestep) of condensation for plumes that have just condensed
               if (cluster.has_condensed() and np.isnan(condensed_time[id])):
                   condensed_time[id] = t
                   cloud = cluster.condensed_mask()
                   #find the minimum height (as a grid coordinate) of condensation
                   min_height = min(cloud)
                   min_height = min_height/(MC['nx']*MC['ny'])
                   min_height = int(min_height)
                   condensed_height[id] = min_height
               #find dry plume region
               #(note that the condensed region is always a subset of the plume in cloudtracker output)
               plume = cluster.plume_mask()
               condensed = cluster.condensed_mask()
               dry_plume = np.setdiff1d(plume, condensed)
               #only increment lifetime of plumes if there are dry plume
               #points (i.e. the plume has not completely formed into a cloud)
              #check here if a cloud has split off from a plume?
               if len(dry_plume) > 0:
                   lifetimes[id] = lifetimes[id] + MC['dt']
                   proj = index_to_xy(dry_plume, MC).T
                   areas[id] = areas[id] + len(proj)
            
    lifetimes = lifetimes[filtered_ids]
    areas = areas[filtered_ids]
    areas = areas*(MC['dx']*MC['dy'])/(lifetimes/MC['dt'])
    #has_condensed = has_condensed[filtered_ids]
    condensed_time = condensed_time[filtered_ids]
    condensed_height = condensed_height[filtered_ids]

    return lifetimes, areas, condensed_time, condensed_height


"""

 inputs: filename  -  name of .pkl file at some time step (contains a dictionary of Cluster objects)
         filtered_ids - ids of clusters after filtering
         MC - a dictionary with the model configuration
         h (optional parameter) - height level
        

 output: an array of distances to plume (dry region only) edge (the index gives the grid point index)

"""

def compute_distances_to_plume_edges(filename, filtered_ids, MC, h=None):
   

    clusters = cPickle.load(open(filename, 'rb'))
    
    plume_depths = label_depths(clusters, filtered_ids, MC, 'plume')
    plume_indices, env_indices = find_plume_indices(clusters, filtered_ids, MC)
    env_halo = find_halo(env_indices, MC)
    env_depth = label_depth(env_indices, env_halo, MC)

    if h == None:
        #label cloud/environment distances to cloud edge
        distances = np.zeros(MC['ny']*MC['nx']*MC['nz'])
        for id, depths in plume_depths.iteritems():
            for n, indices in enumerate(depths):
            #assuming grid spacing is equal in all directions
                distances[indices] = -(n+.5)*MC['dx']
            for n, indices in enumerate(env_depth):
                distances[indices] = (n+.5)*MC['dx']
            return distances
    else:
    #label cloud/environment distances to cloud edge at height level h
        plume_ids_at_h = find_plume_ids_at_z(h, clusters, MC)
        plume_ids_at_h = np.intersect1d(filtered_ids, plume_ids_at_h)
        #if there are clouds at height h, compute distances to cloud edges
        #if len(cloud_ids_at_h):
        plume_depths = label_depths_at_z(h, clusters, plume_ids_at_h, MC, 'plume')
        env_indices = find_indices_at_z(h, env_indices, MC)
        env_halo = find_indices_at_z(h, env_halo, MC)
        env_depth = label_depth(env_indices, env_halo, MC)

        distances_at_h = np.empty(MC['ny']*MC['nx']*MC['nz'])
        distances_at_h[:] = np.NAN


        if len(plume_ids_at_h):
            for id, depths in plume_depths.iteritems():
                for n, indices in enumerate(depths):
                    distances_at_h[indices] = -(n+.5)*MC['dx']

                for n, indices in enumerate(env_depth):
                    distances_at_h[indices] = (n+.5)*MC['dx']

        valid_indices = np.where(np.isfinite(distances_at_h))[0]
        distances_at_h = distances_at_h[valid_indices]

        return distances_at_h, plume_ids_at_h

   








    
    
    







