import numpy as np
import cPickle
import glob
import matplotlib.pyplot as plt


"""
 function to project  an index (from a 3D row-major array) onto the x-y plane
 note: numpy arrays are in row-major order by default
 
 inputs: index - index to project
         ny, nx - grid points in y, and x directions respectively
                      (i.e. the horizontal dimensions of the array)

 output: returns a tuple with the x-y projection of index


"""

def index_to_xy(index, ny, nx):
    index = index % (ny*nx)
    y = index / nx
    x = index % nx
    return (x, y)

"""
 script that outputs plume statistics (mean lifetime, size distribution, etc)

"""

    
filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/clusters*.pkl')
filenames.sort()

dt = 2 #timestep in min
nx = 64 #number of grid points in x and y direction
ny = 64
dx = 50 #grid spacing (m)
dy = 50
clusters_list = []
maxid = 0


for fname in filenames:
    clusters = cPickle.load(open(fname, 'rb'))
    clusters_list.append(clusters)
    ids = clusters.keys()
    maxid_tmp = max(ids)
    if maxid_tmp > maxid:
        maxid = maxid_tmp

#array to keep track of lifetimes of clusters
lifetimes = np.zeros(maxid+1)
#array to keep track of volumes of clusters
#vols = np.zeros(maxid+1)
#array to keep track of areas of xy projections of clusters
areas = np.zeros(maxid+1)



#clusters_list contains the clusters at each time step
#clusters is a dictionary, with cluster id as key,
#and a dictionary as value (containing information about the cluster)
for clusters in clusters_list:
    #iterate over each cluster at the current time
    #increment the lifetime, and find the area of the xy projection
    for id, cluster in clusters.iteritems():
         lifetimes[id] = lifetimes[id] + dt
         for index in cluster['plume']:
             xy_proj = set()
             proj = index_to_xy(index, ny, nx)
             xy_proj.add(proj)
             areas[id] = areas[id] + len(xy_proj)


#filter out noise (i.e. clusters that exist only for a single time step)
plume_ids = np.where(lifetimes > dt)
lifetimes = lifetimes[plume_ids]
areas = areas[plume_ids]
#average the plume projection areas over their lifetime
areas = areas*(dx*dy)/(lifetimes/dt)

print "number of plumes: ",  len(lifetimes)

mean_lt = np.mean(lifetimes)
print "mean plume lifetime (min): ", mean_lt

mean_area = np.mean(areas)
print "mean plume projection area (km^2): ", mean_area/(1e6)

l = np.sqrt(areas)
print "mean length scale (m): ", np.mean(l)

bins = np.arange(31)
plt.figure(1)
plt.hist(lifetimes, bins)
plt.xlabel('lifetime (min)')
plt.ylabel('number of plumes')
plt.figure(2)
bins = np.linspace(0,1000,10)
plt.hist(l, bins)
plt.xlabel(r'projected area$^{1/2}$ (m)')
plt.ylabel('number of plumes')
plt.show()




    
    
    







