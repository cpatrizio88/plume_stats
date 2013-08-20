from util_functions import find_cloud_ids_at_height, find_indices_at_height, index_to_xy
import cPickle
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker/')
from cloudtracker.utility_functions import find_halo

filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/cluster_objects*.pkl')
filenames.sort()

cluster = cPickle.load(open(filenames[0], 'rb'))[0]
MC = cluster.MC

t = 90
h = 720/MC['dz']
clusters = cPickle.load(open(filenames[t], 'rb'))
cloud_edges_at_h = {}
cloud_ids_at_h = find_cloud_ids_at_height(h, clusters, MC)

for id, cluster in clusters.iteritems():
    if id in cloud_ids_at_h:
        halo = find_halo(cluster.condensed_mask(), MC)
        halo_at_h = find_indices_at_height(h, halo, MC)
        halo_at_h = index_to_xy(halo_at_h, MC)
        cloud_edges_at_h[id] = halo_at_h
        
print 'at height {0} m.'.format(h*MC['dz'])
num_clouds = len(cloud_ids_at_h)
print 'number of clouds: ', num_clouds
print 'cloud ids: ', cloud_ids_at_h
colormap = cm.Paired(np.linspace(0, 1, num_clouds))
colors = iter(colormap)

for id, cloud_edge in cloud_edges_at_h.iteritems():
    x = cloud_edge[0, ...]*MC['dx']
    y = cloud_edge[1, ...]*MC['dy']
    plt.scatter(x, y, color = next(colors))
    
plt.title('cloud edges at height {0} m'.format(h*MC['dz']))
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()
