from plume_stats import compute_plume_vars
from cloud_stats import compute_cloud_vars, compute_distances_to_cloud_edges
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util_functions import find_indices_at_z, filter_clusters, find_plume_ids_at_t
from netCDF4 import Dataset
import cPickle
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
from cloudtracker.utility_functions import index_to_zyx


filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/cluster_objects*.pkl')
filenames.sort()

filtered_ids, maxid, MC = filter_clusters(filenames)

t = 85
h = 720/MC['dz']
y = 45

#------plume output-------

lifetimes, areas, condensed_time, condensed_height  = compute_plume_vars(filenames, filtered_ids, maxid, MC)

ismoist = np.isfinite(condensed_time)

condensed_time = condensed_time[ismoist]
condensed_height = condensed_height[ismoist]

lifetimes_dry = lifetimes[~ismoist]/60.
lifetimes_moist = lifetimes[ismoist]/60.

#plume_ids = find_plume_ids_at_t(t, filenames, filtered_ids)

print "number of plumes: ",  len(lifetimes)
print 'number of dry plumes: ', len(lifetimes_dry)
#print 'number of plumes at timestep {0}: {1}'.format(t, len(plume_ids))
print 'number of moist plumes: ', len(lifetimes_moist)
print "dry plume mean lifetime (min): %4.3f " % (np.mean(lifetimes_dry))
dl = np.sqrt(areas[~ismoist])
print "dry plume mean horizontal length scale (m): %4.3f" % (np.mean(dl))
ml = np.sqrt(areas[ismoist])
print "moist plume mean lifetime (min): %4.3f " % (np.mean(lifetimes_moist))
print "moist plume mean horizontal length scale (m): %4.3f" % (np.mean(ml))
print "mean LCL: {0:3.2f} m".format(np.mean(condensed_height)*MC['dz'])

bins = np.arange(31)
plt.hist(lifetimes_dry, bins)
plt.xlabel('lifetime (min)')
plt.ylabel('number of dry plumes')
plt.figure()
bins = np.arange(11)
plt.hist(lifetimes_moist, bins)
plt.xlabel('lifetime (min)')
plt.ylabel('number of moist plumes')
plt.figure()
bins = np.linspace(0,500,50)
plt.hist(ml, bins)
plt.xlabel(r'projected area$^{1/2}$ (m)')
plt.ylabel('number of moist plumes')
plt.figure()
bins = np.linspace(0,500,50)
plt.hist(dl, bins)
plt.xlabel(r'projected area$^{1/2}$ (m)')
plt.ylabel('number of dry plumes')

#------cloud output-------

lifetimes, areas = compute_cloud_vars(filenames, filtered_ids, maxid, MC)

distances = compute_distances_to_cloud_edges(filenames[t], filtered_ids, MC)
distances_at_h, cloud_ids_at_h = compute_distances_to_cloud_edges(filenames[t], filtered_ids, MC, h)


print "number of clouds: ",  len(lifetimes)

print "number of clouds at height {0} m: {1}".format(h*MC['dz'], len(cloud_ids_at_h))

mean_lt = np.mean(lifetimes)
print "cloud mean lifetime (min): %4.3f " % (mean_lt/60.)

l = np.sqrt(areas)
print "cloud mean horizontal length scale (m): %4.3f" % np.mean(l)

#get variable fields at timestep t
fields = glob.glob('/tera/phil/sam_cpatrizi/OUT_3D/BOMEXTR*.nc')
fname = fields[t]
ncfile = Dataset(fname, 'r')
qv = ncfile.variables['QV'][:]
qn = ncfile.variables['QN'][:]
qt = (qn + qv).flatten()


indices_at_h = find_indices_at_z(h, np.arange(len(distances)), MC)

r1, r2 = min(distances), max(distances)
q1, q2 = min(qt), max(qt)
xbins = np.linspace(r1, r2, 101)
ybins = np.linspace(q1, q2, 101)
plt.figure()
hist2d, ybins, xbins = np.histogram2d(qt, distances, bins = (ybins, xbins))
plt.pcolor(xbins, ybins, hist2d, cmap = cm.spectral_r, vmin=0)
plt.colorbar()
plt.xlabel('distance from cloud edge (m)')
plt.ylabel('total mixing ratio (g/kg)')
plt.title('using cloudtracker output at timestep {0}'.format(t))
if len(distances_at_h):
    plt.figure()
    r1, r2 = min(distances_at_h), max(distances_at_h)
    xbins = np.linspace(r1, r2, 101)
    hist2d, ybins, xbins = np.histogram2d(qt[indices_at_h], distances_at_h, bins = (ybins, xbins))
    plt.pcolor(xbins, ybins, hist2d, cmap = cm.spectral_r, vmin=0)
    plt.colorbar()
    plt.xlabel('horizontal distance from cloud edge (m)')
    plt.ylabel('total mixing ratio (g/kg)')
    plt.title('using cloudtracker output at timestep {0}, height {1} m.'.format(t, h*MC['dz']))
bins = np.arange(31)
plt.figure()
plt.hist(lifetimes/60., bins)
plt.xlabel('lifetime (min)')
plt.ylabel('number of clouds')
bins = np.linspace(0,500,50)
plt.figure()
plt.hist(l, bins)
plt.xlabel(r'projected area$^{1/2}$ (m)')
plt.ylabel('number of clouds')
plt.show()
