from util_functions import find_cloud_ids_at_height, find_indices_at_height, index_to_xy, filter_clusters
from cloud_stats import compute_distances_to_cloud_edges
import cPickle
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker/')
from cloudtracker.utility_functions import find_halo


"""
plots cloud edges at height level h

"""

def plot_cloud_edges(filename, filtered_ids, MC,  h):

    clusters = cPickle.load(open(filename, 'rb'))
    cloud_edges_at_h = {}
    cloud_ids_at_h = find_cloud_ids_at_height(h, clusters, MC)
    cloud_ids_at_h = np.intersect1d(cloud_ids_at_h, filtered_ids)

    for id, cluster in clusters.iteritems():
        if id in cloud_ids_at_h:
            halo = find_halo(cluster.condensed_mask(), MC)
            halo_at_h = find_indices_at_height(h, halo, MC)
            halo_at_h = index_to_xy(halo_at_h, MC)
            cloud_edges_at_h[id] = halo_at_h

    num_clouds = len(cloud_ids_at_h)
    colormap = cm.Paired(np.linspace(0, 1, num_clouds))
    colors = iter(colormap)

    for id, cloud_edge in cloud_edges_at_h.iteritems():
        x = cloud_edge[0, ...]*MC['dx']
        y = cloud_edge[1, ...]*MC['dy']
        plt.scatter(x, y, color = next(colors))


"""
plots distance to cloud edge (contour plot) at height level h

"""

def plot_cloud_edge_distance(filename, filtered_ids, MC, h):

    nx = MC['nx']
    ny = MC['ny']

    x = np.arange(nx)
    x = x*MC['dx']
    y = np.arange(ny)
    y = y*MC['dy']
    x, y = np.meshgrid(x, y)

    distances_at_h, cloud_ids_at_h = compute_distances_to_cloud_edges(filename, filtered_ids, MC, h)

    levs = np.linspace(np.min(distances_at_h), np.max(distances_at_h), 20)
    distances_at_h = distances_at_h.reshape(nx, ny)
    plt.contourf(x, y, distances_at_h, levels=levs)
    plt.colorbar()  

    
