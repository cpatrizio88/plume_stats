from mpl_toolkits.mplot3d import Axes3D
from plume_stats import plume_stats
from cloud_stats import cloud_stats
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util_functions import find_indices_at_height
from netCDF4 import Dataset
import cPickle
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker

#------plume output-------

filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/cluster_objects*.pkl')
filenames.sort()
cluster = cPickle.load(open(filenames[0], 'rb'))[0]
#model configuration (grid spacing, domain size, etc)
MC = cluster.MC

lifetimes, areas, has_condensed = plume_stats(filenames)


#dry_plume_ids = np.where(~has_condensed)

print "number of plumes: ",  len(lifetimes)
print 'number of dry plumes: ', len(lifetimes[~has_condensed])
print 'number of moist plumes: ', len(lifetimes[has_condensed])

#mean_lt = np.mean(lifetimes[plume_ids])
#print "mean plume lifetime (min): %4.3f " % (mean_lt)
print "moist plume mean lifetime (min): %4.3f " % (np.mean(lifetimes[has_condensed]))
print "dry plume mean lifetime (min): %4.3f " % (np.mean(lifetimes[~has_condensed]))

ml = np.sqrt(areas[has_condensed])
dl = np.sqrt(areas[~has_condensed])
print "moist plume mean horizontal length scale (m): %4.3f" % (np.mean(ml))
print "dry plume mean horizontal length scale (m): %4.3f" % (np.mean(dl))

bins = np.arange(31)
plt.hist(lifetimes, bins)
plt.xlabel('lifetime (min)')
plt.ylabel('number of plumes')
plt.figure()
bins = np.arange(11)
plt.hist(lifetimes[has_condensed], bins)
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
t = 90
lifetimes, areas, distances = cloud_stats(filenames, t)

print "number of clouds: ",  len(lifetimes)

mean_lt = np.mean(lifetimes)
print "cloud mean lifetime (min): %4.3f " % (mean_lt)

l = np.sqrt(areas)
print "cloud mean horizontal length scale (m): %4.3f" % np.mean(l)

#get variable fields at timestep t
t = 90
fields = glob.glob('/tera/phil/sam_cpatrizi/OUT_3D/BOMEXTR*.nc')
fname = fields[t]
ncfile = Dataset(fname, 'r')
qv = ncfile.variables['QV'][:]
qn = ncfile.variables['QN'][:]
qt = (qn + qv).flatten()

#look at a single height level
h = 720/MC['dz']

indices_at_h = find_indices_at_height(h, np.arange(len(distances)), MC)
distances_at_h = distances[indices_at_h]

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
plt.figure()
r1, r2 = min(distances_at_h), max(distances_at_h)
xbins = np.linspace(r1, r2, 101)
hist2d, ybins, xbins = np.histogram2d(qt[indices_at_h], distances_at_h, bins = (ybins, xbins))
plt.pcolor(xbins, ybins, hist2d, cmap = cm.spectral_r, vmin=0)
plt.colorbar()
plt.xlabel('distance from cloud edge (m)')
plt.ylabel('total mixing ratio (g/kg)')
plt.title('using cloudtracker output at timestep {0}, height {1} m.'.format(t, h*MC['dz']))

bins = np.arange(31)
plt.figure()
plt.hist(lifetimes, bins)
plt.xlabel('lifetime (min)')
plt.ylabel('number of clouds')
bins = np.linspace(0,500,50)
plt.hist(l, bins)
plt.xlabel(r'projected area$^{1/2}$ (m)')
plt.ylabel('number of clouds')
plt.show()
