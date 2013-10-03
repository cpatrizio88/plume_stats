import glob
import cPickle
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from plot_clusters import contour_cloud_edge_distance, plot_cloud_edges_z, plot_plumes_xz, \
                          plot_plumes_y, plot_plumes_xy, plot_plumes_z, contour_plume_edge_distance 
from util_functions import filter_clusters, find_cloud_ids_at_z, find_plume_ids_at_z

filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/cluster_objects*.pkl')
filenames.sort()

filtered_ids, maxid, MC = filter_clusters(filenames)

t = 68
h = 720/MC['dz']

clusters = cPickle.load(open(filenames[t], 'rb'))

for z in range(MC['nz']):
   plume_ids = find_plume_ids_at_z(z, clusters, MC)
   cloud_ids = find_cloud_ids_at_z(z, clusters, MC)
   #print cloud_ids
   if len(plume_ids) or len(cloud_ids):
      print 'height level: {0}, number of plumes: {1}, number of clouds: {2}'.format(z, len(plume_ids), len(cloud_ids))


z = 11
plt.figure()
contour_plume_edge_distance(filenames[t], filtered_ids, MC, z)
plot_plumes_z(filenames[t], filtered_ids, MC, z)
plt.title('plume (dry region only) horizontal cross section at height {0} m, timestep {1}'.format(z*MC['dz'], t))
plt.figure()
plot_plumes_xz(filenames[t], filtered_ids, MC, 'add clouds')
plt.title('x-z projections of plumes (w/ clouds) at timestep {0}'.format(t))
#plt.figure()
#plot_plumes_xz(filenames[t], filtered_ids, MC)
#plt.title('x-z projections of plumes (dry region only) at timestep {0}'.format(t))
#plt.figure()
#plot_plumes_y(filenames[t], filtered_ids, MC, y)
#plt.ylim((0,1500))
#plt.xlim((0, MC['nx']*MC['dx']))
#plt.figure()
#plot_plumes_xy(filenames[t], filtered_ids, MC)
#plt.title('x-y projections of plumes (dry region only) at timestep {0}'.format(t))
plt.figure()
contour_cloud_edge_distance(filenames[t], filtered_ids, MC, h)
plot_cloud_edges_z(filenames[t], filtered_ids, MC, h)
plt.title('distance to cloud edges at height {0} m, timestep {1}'.format(h*MC['dz'], t))
plt.show()
