import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
import cPickle
import glob
from netCDF4 import Dataset
from util_functions import index_to_xy


"""
  outputs statistics for plumes stored as Cluster objects in .pkl files (cloudtracker output)

  input: filenames - contains list names of .pkl files 

  outputs: array of mean lifetime (min)
           arary of mean size (m)
           array of booleans indicating whether a plume forms a condensed region or not

           note: all plumes that exist for only a single model output time step are ignored 
"""

def plume_stats(filenames):

    filenames.sort()
    cluster = cPickle.load(open(filenames[0], 'rb'))[0]
    MC = cluster.MC

    #get grid parameters
    dt = MC['dt']/60. #timestep in min
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
    #lifetimes of clusters before condensing
    lt_before_condensed = np.zeros(maxid+1)
    
    #areas of xy projections of plumes (before condensed region forms)
    areas = np.zeros(maxid+1)
    #keep track of plumes that form condensed regions
    #initialize all plumes to False
    has_condensed = np.zeros(maxid+1) < -1
    
    #clusters is a dictionary of Cluster objects
    #(see cloud_objects.py in cloudtracker)
    for t, clusters in enumerate(clusters_list):
        #iterate over each cluster at the current time
        #increment the lifetime of the plume
        #and find the area of plume x-y projection
        for id, cluster in clusters.iteritems():
           lifetimes[id]  = lifetimes[id] + dt
           if not has_condensed[id]:
             lt_before_condensed[id] = lt_before_condensed[id] + dt
             proj = index_to_xy(cluster.plume_mask(), MC).T
             areas[id] = areas[id] + len(proj)
           #if plume has condensed
           if cluster.has_condensed():
               has_condensed[id] = True

    #filter out noise (i.e. clusters that exist only for a single time step)
    plume_ids = np.where(lifetimes > dt)[0]

    #average the plume projection areas over their lifetime
    areas = areas*(dx*dy)/(lt_before_condensed/dt)

    lifetimes = lifetimes[plume_ids]
    lt_before_condensed = lt_before_condensed[plume_ids]
    areas = areas[plume_ids]
    has_condensed = has_condensed[plume_ids]

    return lt_before_condensed, areas, has_condensed

   








    
    
    







