from __future__ import division
import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
import cPickle
import glob
from netCDF4 import Dataset
from util_functions import index_to_xy, filter_clusters


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
    #a NAN entry signifies that the plume has not condensed
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




   








    
    
    







