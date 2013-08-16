import glob
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cPickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import cloudtracker
from cloudtracker.utility_functions import find_halo
from util_functions import xy_to_index, index_to_xy, expand_indexes_horizontal


#script to plot x-y projections of clouds

filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/cluster_objects*.pkl')
filenames.sort()

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

#keep track of x-y projections of clouds at all times
projections = []

MC = cPickle.load(open(filenames[0], 'rb'))[0].MC

#clusters is a dictionary of Cluster objects
#(see cloud_objects.py in cloudtracker)
for clusters in clusters_list:
   #keep track of x-y projections of clouds at current time step
   cloudprojs = {}
   for id, cluster in clusters.iteritems():
      if cluster.has_condensed():
         lifetimes[id] = lifetimes[id] + MC['dt']
         proj = index_to_xy(cluster.condensed_mask(), MC)
         cloudprojs[id] = proj
   projections.append(cloudprojs)

cloud_ids = np.where(lifetimes > MC['dt'])[0]

t = 90
projs = projections[t]
ids = projs.keys()
cloud_ids = np.intersect1d(ids, cloud_ids)

num_clouds = len(cloud_ids)

print 'number of clouds: ', num_clouds
print 'cloud ids: ', cloud_ids
colormap = cm.Paired(np.linspace(0, 1, num_clouds))
colors = iter(colormap)

ax1 = plt.gca()

for id, cloud in projs.iteritems():
    if id in cloud_ids:
        x = cloud[0, ...]
        y = cloud[1, ...]
        indices = xy_to_index(x, y, MC)
        expanded_indices = expand_indexes_horizontal(indices, MC)
        boundary = np.setdiff1d(expanded_indices, indices)
        boundary = index_to_xy(boundary, MC)
        x = boundary[0, ...]
        y = boundary[1, ...]
        plt.scatter(x, y, color = next(colors), s = 7*np.ones(x.shape))
        
ax1.set_xlim((-5,70))
ax1.set_ylim((-5,70))
plt.show()      

            






    
        
        
