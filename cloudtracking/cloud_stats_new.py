import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
from cloudtracker.utility_functions import index_to_zyx, zyx_to_index, find_halo, expand_indexes
import cPickle
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from netCDF4 import Dataset
from pylab import *

"""
 script that outputs cloud statistics (mean lifetime, size distribution, etc)

"""

def main():
    
    filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/cluster_objects*.pkl')
    filenames.sort()

    cluster = cPickle.load(open(filenames[0], 'rb'))[0]

    MC = cluster.MC

    #get grid parameters
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
    #areas of xy projections of clusters
    areas = np.zeros(maxid+1)
    
    #clusters is a dictionary of Cluster objects
    #(see cloud_objects.py in cloudtracker)
    for t, clusters in enumerate(clusters_list):
        #iterate over each cluster at the current time
        #increment the lifetime if it has a condensed region,
        # and find the area of xy cloud projection
        for id, cluster in clusters.iteritems():
            if cluster.has_condensed():
                xy_proj = set()
                lifetimes[id] = lifetimes[id] + dt
            for index in cluster.condensed_mask():
                proj = index_to_xy(index, ny, nx)
                xy_proj.add(proj)
            areas[id] = areas[id] + len(xy_proj)

    #filter out noise (i.e. clusters that exist only for a single time step)
    cloud_ids = np.where(lifetimes > dt)[0]
    
    t = 90
    clusters = clusters_list[t]
    #find cloud/environment depths for each cloud at time step t
    cloud_depths = find_cloud_depths(clusters, cloud_ids, MC)
    env_indices, cloud_indices = find_all_indices(clusters, MC)
    env_halo = find_halo(env_indices, MC)
    env_depth = label_depth(env_indices, env_halo, MC)

    lifetimes = lifetimes[cloud_ids]
    areas = areas[cloud_ids]
    #average the cloud projection areas over their lifetime
    areas = areas*(dx*dy)/(lifetimes/dt)

    print "number of clouds: ",  len(lifetimes)

    mean_lt = np.mean(lifetimes)
    print "mean cloud lifetime (min): %4.3f " % (mean_lt)

    mean_area = np.mean(areas)
    #print "mean cloud projection area (km^2): %4.3f " % (mean_area/(1e6))

    l = np.sqrt(areas)
    print "mean horizontal length scale (m): %4.3f" % (np.mean(l))

    print "cloud ids at time step %2d: " % (t), cloud_depths.keys()

    #test at a single height level
    h = 760/dz

    #get variable fields at timestep t
    filenames = glob.glob('/tera/phil/sam_cpatrizi/OUT_3D/BOMEXTR*.nc')
    fname = filenames[t]
    ncfile = Dataset(fname, 'r')
    qv = ncfile.variables['QV'][:]
    qn = ncfile.variables['QN'][:]
    qt = (qn + qv).flatten()

    #distances contains grid point distance (in m) from cloud boundary
    distances = np.zeros_like(qt)
     
    print 'calculating cloud grid point distances to cloud edges'
    for id, depths in cloud_depths.iteritems():
       for n, indices in enumerate(depths):
            #assuming grid spacing is equal in all directions
            distances[indices] = -(n+.5)*MC['dx']
            
    print 'calculating environment grid point distances to cloud edges'
    for n, indices in enumerate(env_depth):
           distances[indices] = (n+.5)*MC['dx']
    
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
    plt.figure()
    r1, r2 = min(distances_at_h), max(distances_at_h)
    xbins = np.linspace(r1, r2, 101)
    hist2d, ybins, xbins = np.histogram2d(qt[indices_at_h], distances_at_h, bins = (ybins, xbins))
    plt.pcolor(xbins, ybins, hist2d, cmap = cm.spectral_r, vmin=0)
    plt.colorbar()
    plt.xlabel('distance from cloud edge (m)')
    plt.ylabel('total mixing ratio (g/kg)')
    plt.title('at height {0} m.'.format(h*dz))
    
    bins = np.arange(31)
    plt.figure()
    plt.hist(lifetimes, bins)
    plt.xlabel('lifetime (min)')
    plt.ylabel('number of clouds')
    plt.figure()
    bins = np.linspace(0,1000,50)
    plt.hist(l, bins)
    plt.xlabel(r'projected area$^{1/2}$ (m)')
    plt.ylabel('number of clouds')
    plt.show()

"""
 function to project an index (from a 3D row-major array) onto the x-y plane
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

returns a list of arrays, each array contains the grid indices for a depth level
in a cloud/cloud environment 

based off of label_distance bomex_histogram_data.py  in /users/jdawe/figures/knut

inputs: indices - grid indices specifying the cloud/environment
        MC - dictionary containing grid parameters (specified in model_config.cfg)
        halo - grid indices adjacent to the cloud boundary

output: a list of numpy arrays with grid indices at each depth level (grid indices
        that are equidistant from the cloud boundary) in a cloud/cloud environment
        (starts at grid indices closest to the cloud boundary) 
        
"""

def label_depth(indices, halo, MC):

   mask_indexes = indices[:]

   done = False
   depth = []
   while not done:
        expanded_indexes = expand_indexes(halo, MC)

        # From the expanded indexes list, select the points in the indexes list
        # expand_index_list returns unique values so we don't have to filter duplicates
        new_indexes = np.intersect1d(expanded_indexes, mask_indexes, assume_unique=True)

        if len(new_indexes) == 0:
            done = True
        else:
            depth.append(new_indexes)

            # From the mask_indexes, remove all the new points
            mask_indexes = np.setdiff1d(mask_indexes, new_indexes, assume_unique=True)

            halo = new_indexes

   # return a list composed of arrays of the points at each depth
   return depth


"""
 returns a dictionary (indexed by cloud id) containing cloud depth arrays
 (as described in label_depth) 

inputs: clusters - a dictionary of Cluster objects at a single time step (indexed by cloud id)
         cloud_ids - a list of cloud ids, must be a subset of the id's found in clusters
         MC - dictionary containing grid parameters (as in model_config.cfg)

 outputs: see above

"""

def find_cloud_depths(clusters, cloud_ids, MC):
    cloud_depths = {}
    env_indices, cloud_indices = find_all_indices(clusters, MC)
    for id, cluster in clusters.iteritems():
        if id in cloud_ids:
              indices = cluster.condensed_mask()
              cloud_halo = cluster.condensed_halo()
              cloud_depths[id] = label_depth(indices, cloud_halo, MC)

    return cloud_depths
              

"""
  returns all grid indices that are not condensed points
  and all grid indices that are condensed points as numpy arrays

"""

def find_all_indices(clusters, MC):
    
    cloud_indices = np.array([])

    grid = np.arange(MC['nx']*MC['ny']*MC['nz'])

    for id, cluster in clusters.iteritems():
        tmp = cluster.condensed_mask()
        cloud_indices = np.concatenate([cloud_indices, tmp])

    env_indices = np.setdiff1d(grid, cloud_indices, assume_unique=True)

    return env_indices, cloud_indices

#calculates the distance between index1 and index2
#(note: index2 or index1 may be an array of indices, but not both)
def calc_distance(index1, index2, MC):

    point1 = index_to_zyx(index1, MC)
    
    point2 = index_to_zyx(index2, MC)
    
    # Calculate distances
    ny, nx = MC['ny'], MC['nx']
    dy, dx, dz = MC['dy'], MC['dx'], MC['dz']
        
    delta_x = (point2[2, ...] - point1[2, ...])*dx

    delta_y = (point2[1, ...] - point1[1, ...])*dy
   
    delta_z =(point2[0, ...] - point1[0, ...])*dz

    return np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)

def find_indices_at_height(h, indices, MC):

    nx = MC['nx']
    ny = MC['ny']
    
    hit = np.logical_and(h*nx*ny < indices, indices <= (h+1)*nx*ny)
    return indices[hit]

if __name__ == '__main__':
    main()






    
    
    







