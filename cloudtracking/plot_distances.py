from cloud_stats import cloud_stats
import matplotlib.pyplot as plt
import glob
import cPickle
from util_functions import index_to_xy
import numpy as np

filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/cluster_objects*.pkl')
filenames.sort()

cluster = cPickle.load(open(filenames[0], 'rb'))[0]
MC = cluster.MC
t=90
h=720/MC['dz']
nx = MC['nx']
ny = MC['ny']

x = np.arange(nx)
x = x*MC['dx']
y = np.arange(ny)
y = y*MC['dy']
x, y = np.meshgrid(x, y)

lifetimes, areas, distances, distances_at_h, cloud_ids_at_h = cloud_stats(filenames, t, h)

print 'number of clouds', len(cloud_ids_at_h)

levs = np.linspace(np.min(distances_at_h) ,np. max(distances_at_h), 20)
distances_at_h = distances_at_h.reshape(nx, ny)
c = plt.contourf(x, y, distances_at_h, levels=levs)
plt.colorbar()
plt.show()
plt.title('distance to cloud edge (m) at height {0} m, timestep {1}'.format(h*MC['dz'], t))
plt.xlabel('x (m)')
plt.ylabel('y (m)')

