import cPickle
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker/')
from cloudtracker.utility_functions import find_halo
from cloud_stats import compute_distances_to_cloud_edges
from plume_stats import compute_distances_to_plume_edges
from util_functions import find_cloud_ids_at_z, find_indices_at_z, \
                           index_to_xy, filter_clusters, index_to_xz, find_plume_ids_at_y, \
                           find_indices_at_y, find_plume_ids_at_z


"""

plots a horizontal slice of cloud edges at height level h

"""
def plot_cloud_edges_z(filename, filtered_ids, MC,  h):

    clusters = cPickle.load(open(filename, 'rb'))
    cloud_edges = {}
    
    cloud_ids = find_cloud_ids_at_z(h, clusters, MC)
    cloud_ids = np.intersect1d(cloud_ids, filtered_ids)

    for id, cluster in clusters.iteritems():
        if id in cloud_ids:
            halo = find_halo(cluster.condensed_mask(), MC)
            halo_at_h = find_indices_at_z(h, halo, MC)
            halo_at_h = index_to_xy(halo_at_h, MC)
            cloud_edges[id] = halo_at_h

    num_clouds = len(cloud_ids)
    colormap = cm.Paired(np.linspace(0, 1, num_clouds))
    colors = iter(colormap)

    for id, cloud_edge in cloud_edges.iteritems():
        x = cloud_edge[0, ...]*MC['dx']
        y = cloud_edge[1, ...]*MC['dy']
        plt.scatter(x, y, color = next(colors), s = 5*np.ones(x.shape))

    ax1 = plt.gca()
    ax1.set_xlim(0, MC['nx']*MC['dx'])
    ax1.set_ylim(0, MC['ny']*MC['dy'])



def plot_plumes_z(filename, filtered_ids, MC, h):

    clusters = cPickle.load(open(filename, 'rb'))
    plumes = {}

    plume_ids = find_plume_ids_at_z(h, clusters, MC)
    plume_ids = np.intersect1d(plume_ids, filtered_ids)

    for id, cluster in clusters.iteritems():
        if id in plume_ids:
           dry_plume = np.setdiff1d(cluster.plume_mask(), cluster.condensed_mask())
           dry_plume = find_indices_at_z(h, dry_plume, MC)
           if len(dry_plume):
               dry_plume = index_to_xy(dry_plume, MC)
               plumes[id] = dry_plume

    num_plumes = len(plume_ids)
    colormap = cm.Paired(np.linspace(0, 1, num_plumes))
    colors = iter(colormap)

    for id, plume in plumes.iteritems():
        x = plume[0, ...]*MC['dx']
        y = plume[1, ...]*MC['dy']
        plt.scatter(x, y, color = next(colors), s = 5*np.ones(x.shape))

    ax1 = plt.gca()
    ax1.set_xlim(0, MC['nx']*MC['dx'])
    ax1.set_ylim(0, MC['ny']*MC['dy'])  
   
    
"""

plots x-z projections of plumes 

if optional parameter is specified (any value):  plot entire plume (including cloud)

if optional parameter not specified: plot plume dry region only

"""
def plot_plumes_xz(filename, filtered_ids, MC, *args):
    if len(args):
        clusters = cPickle.load(open(filename, 'rb'))
        plume_projs = {}

        for id, cluster in clusters.iteritems():
            if id in filtered_ids:
                #halo = find_halo(cluster.plume_mask(), MC)
                plume = cluster.plume_mask()
                if len(plume):
                   plume_proj = index_to_xz(plume, MC)
                   plume_projs[id] = plume_proj

        num_plumes = len(plume_projs.keys())

        colormap = cm.Paired(np.linspace(0, 1, num_plumes))
        colors = iter(colormap)

        for id, plume_proj in plume_projs.iteritems():
            x = plume_proj[0, ...]*MC['dx']
            z = plume_proj[1, ...]*MC['dz']
            plt.scatter(x, z, color = next(colors), s = 8*np.ones(x.shape))

        ax1 = plt.gca()    
        ax1.set_ylim((0,MC['nz']*MC['dz']/2.))
        ax1.set_xlim((0, MC['nx']*MC['dx']))

    else:
        clusters = cPickle.load(open(filename, 'rb'))
        plume_projs = {}

        for id, cluster in clusters.iteritems():
            if id in filtered_ids:
                #halo = find_halo(cluster.plume_mask(), MC)
                dry_plume = np.setdiff1d(cluster.plume_mask(), cluster.condensed_mask())
                if len(dry_plume):
                   plume_proj = index_to_xz(dry_plume, MC)
                   plume_projs[id] = plume_proj

        num_plumes = len(plume_projs.keys())

        colormap = cm.Paired(np.linspace(0, 1, num_plumes))
        colors = iter(colormap)

        for id, plume_proj in plume_projs.iteritems():
            x = plume_proj[0, ...]*MC['dx']
            z = plume_proj[1, ...]*MC['dz']
            plt.scatter(x, z, color = next(colors), s = 8*np.ones(x.shape))

        ax1 = plt.gca()    
        ax1.set_ylim((0,MC['nz']*MC['dz']/2.))
        ax1.set_xlim((0, MC['nx']*MC['dx']))

"""

plots x-y projections of plumes (dry region only)

"""
def plot_plumes_xy(filename, filtered_ids, MC):

    clusters = cPickle.load(open(filename, 'rb'))
    plume_projs = {}

    for id, cluster in clusters.iteritems():
        if id in filtered_ids:
            #halo = find_halo(cluster.plume_mask(), MC)
            dry_plume = np.setdiff1d(cluster.plume_mask(), cluster.condensed_mask())
            if len(dry_plume):
               plume_proj = index_to_xy(dry_plume, MC)
               plume_projs[id] = plume_proj
 
    num_plumes = len(plume_projs.keys())

    colormap = cm.Paired(np.linspace(0, 1, num_plumes))
    colors = iter(colormap)

    for id, plume_proj in plume_projs.iteritems():
        x = plume_proj[0, ...]*MC['dx']
        y = plume_proj[1, ...]*MC['dy']
        plt.scatter(x, y, color = next(colors), s = 8*np.ones(x.shape))

    ax1 = plt.gca()
    ax1.set_xlim(0, MC['nx']*MC['dx'])
    ax1.set_ylim(0, MC['ny']*MC['dy'])
    
"""

plots plumes (dry region only) at horizontal level y

"""
def plot_plumes_y(filename, filtered_ids, MC, y):

    clusters = cPickle.load(open(filename, 'rb'))
    plumes = {}

    plume_ids = find_plume_ids_at_y(y, clusters, MC)
    plume_ids = np.intersect1d(filtered_ids, plume_ids)

    for id, cluster in clusters.iteritems():
        if id in plume_ids:
            dry_plume = np.setdiff1d(cluster.plume_mask(), cluster.condensed_mask())
            dry_plume = find_indices_at_y(y, dry_plume, MC)
            dry_plume = index_to_xz(dry_plume, MC)
            plumes[id] = dry_plume

    num_plumes = len(plumes.keys())

    colormap = cm.Paired(np.linspace(0, 1, num_plumes))
    colors = iter(colormap)

    for id, plume in plumes.iteritems():
        x = plume[0, ...]*MC['dx']
        z = plume[1, ...]*MC['dz']
        plt.scatter(x, z, color = next(colors), s = 8*np.ones(x.shape))


"""
contour distance to cloud edge (contour plot) at height level h

"""
def contour_cloud_edge_distance(filename, filtered_ids, MC, h):

    nx = MC['nx']
    ny = MC['ny']

    x = np.arange(nx)
    x = x*MC['dx']
    y = np.arange(ny)
    y = y*MC['dy']
    x, y = np.meshgrid(x, y)

    distances_at_h, cloud_ids_at_h = compute_distances_to_cloud_edges(filename, filtered_ids, MC, h)
    print len(distances_at_h)
    levs = np.linspace(np.min(distances_at_h), np.max(distances_at_h), 40)

    
    distances_at_h = distances_at_h.reshape(nx, ny)
    
    plt.contourf(x, y, distances_at_h, levels=levs)
    plt.colorbar()

    ax1 = plt.gca()
    ax1.set_xlim(0, MC['nx']*MC['dx'])
    ax1.set_ylim(0, MC['ny']*MC['dy'])

"""
contour distance to cloud edge (contour plot) at height level h

"""
def contour_plume_edge_distance(filename, filtered_ids, MC, h):

    nx = MC['nx']
    ny = MC['ny']

    x = np.arange(nx)
    x = x*MC['dx']
    y = np.arange(ny)
    y = y*MC['dy']
    x, y = np.meshgrid(x, y)

    distances_at_h, plume_ids_at_h = compute_distances_to_plume_edges(filename, filtered_ids, MC, h)\

    levs = np.linspace(np.min(distances_at_h), np.max(distances_at_h), 40)
    distances_at_h = distances_at_h.reshape(nx, ny)
    
    plt.contourf(x, y, distances_at_h, levels=levs)
    plt.colorbar()

    ax1 = plt.gca()
    ax1.set_xlim(0, MC['nx']*MC['dx'])
    ax1.set_ylim(0, MC['ny']*MC['dy'])


    
