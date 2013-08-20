import numpy as np
import site
import sys
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
from cloudtracker.utility_functions import index_to_zyx, zyx_to_index, find_halo, expand_indexes
import cPickle
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from netCDF4 import Dataset
from pylab import *
from util_functions import index_to_xy, label_cloud_depths, find_all_indices, label_depth, \
                           find_cloud_ids_at_height, find_indices_at_height, \
                           label_cloud_depths_at_height

"""
 outputs cloud statistics (mean lifetime, size, distances to cloud edge)

"""

def cloud_stats(filenames, t, h):
    
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
    #lifetimes of clouds
    cloud_lifetimes = np.zeros(maxid+1)
    #areas of xy projections of clouds
    areas = np.zeros(maxid+1)
    
    #clusters is a dictionary of Cluster objects
    #(see cloud_objects.py in cloudtracker)
    for t, clusters in enumerate(clusters_list):
        #iterate over each cluster at the current time
        #increment the lifetime if it has a condensed region,
        #and find the area of xy cloud projection
        for id, cluster in clusters.iteritems():
            lifetimes[id] = lifetimes[id] + dt
            if cluster.has_condensed():
                cloud_lifetimes[id] = cloud_lifetimes[id] + dt
                proj = index_to_xy(cluster.condensed_mask(), MC).T
                areas[id] = areas[id] + len(proj)

    #filter out noise (i.e. clusters that exist only for a single time step)
    cloud_ids = np.where(np.logical_and(lifetimes > dt, cloud_lifetimes > 0))[0]
    
    clusters = clusters_list[t]
    #find cloud/environment depths at timestep t
    cloud_depths = label_cloud_depths(clusters, cloud_ids, MC)
    env_indices, cloud_indices = find_all_indices(clusters, MC)
    env_halo = find_halo(env_indices, MC)
    env_depth = label_depth(env_indices, env_halo, MC)
    #find cloud/environment depths at height h
    cloud_ids_at_h = find_cloud_ids_at_height(h, clusters, MC)
    if len(cloud_ids_at_h):
        cloud_depths_at_h = label_cloud_depths_at_height(h, clusters, MC)
        env_indices_at_h = find_indices_at_height(h, env_indices, MC)
        env_halo_at_h = find_indices_at_height(h, env_halo, MC)
        env_depth_at_h = label_depth(env_indices_at_h, env_halo_at_h, MC)
        #cloud_indices_at_h = find_indices_at_height(h, cloud_indices, MC)
        #print 'length env_indices_at_h:', len(env_indices_at_h)
        #print 'length cloud_indices_at_h:', len(cloud_indices_at_h)
        #count = 0
        #for depth in env_depth_at_h:
        #    count = count + len(depth)
        #print 'length env_depth_at_h:', count
        #count = 0
        #for id, depths in cloud_depths_at_h.iteritems():
        #    for a in depths:
        #        count = count + len(a)
        #print 'length cloud_depths_at_h:', count
    
    cloud_lifetimes = cloud_lifetimes[cloud_ids]
    areas = areas[cloud_ids]
    #average the cloud projection areas over their lifetime
    areas = areas*(dx*dy)/(cloud_lifetimes/dt)

    distances = np.zeros(ny*nx*nz)
    distances_at_h = np.empty(distances.shape)
    distances_at_h[:] = np.NAN

    for id, depths in cloud_depths.iteritems():
        for n, indices in enumerate(depths):
        #assuming grid spacing is equal in all directions
            distances[indices] = -(n+.5)*MC['dx']

    for n, indices in enumerate(env_depth):
        distances[indices] = (n+.5)*MC['dx']

    if len(cloud_ids_at_h):
        for id, depths in cloud_depths_at_h.iteritems():
            for n, indices in enumerate(depths):
                distances_at_h[indices] = -(n+.5)*MC['dx']

            for n, indices in enumerate(env_depth_at_h):
                distances_at_h[indices] = (n+.5)*MC['dx']

    valid_indices = np.where(np.isfinite(distances_at_h))[0]
    #print valid_indices
    distances_at_h = distances_at_h[valid_indices]

    return cloud_lifetimes, areas, distances, distances_at_h, cloud_ids_at_h





    
    
    







