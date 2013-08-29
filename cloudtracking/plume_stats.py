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

def plume_stats(filenames, filtered_ids, maxid, MC):

    filenames.sort()

    #lifetimes of clusters
    lifetimes = np.zeros(maxid+1)
    #lifetimes of clusters before condensing
    lifetimes_moist = np.zeros(maxid+1)
    
    #areas of xy projections of plumes (before condensed region forms)
    areas = np.zeros(maxid+1)
    #keep track of plumes that form condensed regions
    #initialize all plumes to False
    has_condensed = np.zeros(maxid+1) < -1

    
    for fname in filenames:
        clusters = cPickle.load(open(fname, 'rb'))
        #iterate over each cluster at the current time
        #increment the lifetime of the plume
        #and find the area of plume x-y projection
        for id, cluster in clusters.iteritems():
           if id in filtered_ids:
                #only increment lifetime of plumes if there are dry plume
                #points (i.e. the plume is has not completely formed into a cloud)
                if has_condensed[id]:
                    plume = cluster.plume_mask()
                    #the condensed region is always a subset of the plume 
                    condensed = cluster.condensed_mask()
                    dry_plume = np.setdiff1d(plume, condensed)
                    if len(dry_plume) > 0:
                        lifetimes[id] = lifetimes[id] + MC['dt']
                        proj = index_to_xy(dry_plume, MC).T
                        areas[id] = areas[id] + len(proj)
                else:
                    lifetimes[id] = lifetimes[id] + MC['dt']
                    proj = index_to_xy(cluster.plume_mask(), MC).T
                    areas[id] = areas[id] + len(proj)
           #label plumes that have just condensed, find height and time of condensation
           if cluster.has_condensed() and not has_condensed[id]:
               has_condensed[id] = True

    lifetimes = lifetimes[filtered_ids]
    areas = areas[filtered_ids]
    areas = areas*(MC['dx']*MC['dy'])/(lifetimes/MC['dt'])
    has_condensed = has_condensed[filtered_ids]

    return lifetimes, areas, has_condensed

   








    
    
    







