import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
from cloudtracker.utility_functions import index_to_zyx, zyx_to_index, find_halo, expand_indexes
import cPickle
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from netCDF4 import Dataset
from pylab import *
from util_functions import index_to_xy, label_cloud_depths, find_all_indices, label_depth

"""
 outputs cloud statistics (mean lifetime, size)

"""

def cloud_stats(filenames, t):
    
    filenames.sort()
    cluster = cPickle.load(open(filenames[0], 'rb'))[0]
    MC = cluster.MC
    #get model parameters
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
    #dictionary to keep track of cloud depth (indexed by cloud id)
    cloud_depths = {}
    env_depth = []
    
    for fname in filenames:
        clusters = cPickle.load(open(fname, 'rb'))
        clusters_list.append(clusters)
        ids = clusters.keys()
        maxid_tmp = max(ids)
        if maxid_tmp > maxid:
            maxid = maxid_tmp

    #lifetimes of clusters
    lifetimes = np.zeros(maxid+1)
    #areas of xy projections of clusters
    areas = np.zeros(maxid+1)
    
    #clusters is a dictionary of Cluster objects
    #(see cloud_objects.py in cloudtracker)
    for t, clusters in enumerate(clusters_list):
        #iterate over each cluster at the current time
        #increment the lifetime if it has a condensed region,
        #and find the area of xy cloud projection
        for id, cluster in clusters.iteritems():
            if cluster.has_condensed():
                lifetimes[id] = lifetimes[id] + dt
                proj = index_to_xy(cluster.condensed_mask(), MC).T
                areas[id] = areas[id] + len(proj)

    #filter out noise (i.e. clusters that exist only for a single time step)
    cloud_ids = np.where(lifetimes > dt)[0]
    
    clusters = clusters_list[t]
    #find cloud/environment depths for each cloud at time step t
    cloud_depths = label_cloud_depths(clusters, cloud_ids, MC)
    env_indices, cloud_indices = find_all_indices(clusters, MC)
    env_halo = find_halo(env_indices, MC)
    env_depth = label_depth(env_indices, env_halo, MC)

    lifetimes = lifetimes[cloud_ids]
    areas = areas[cloud_ids]
    #average the cloud projection areas over their lifetime
    areas = areas*(dx*dy)/(lifetimes/dt)

    distances = np.zeros(ny*nx*nz)

    for id, depths in cloud_depths.iteritems():
        for n, indices in enumerate(depths):
        #assuming grid spacing is equal in all directions
            distances[indices] = -(n+.5)*MC['dx']

    for n, indices in enumerate(env_depth):
        distances[indices] = (n+.5)*MC['dx']

    return lifetimes, areas, distances





    
    
    







