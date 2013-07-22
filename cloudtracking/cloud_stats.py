import numpy as np
import site
site.addsitedir('/home/cpatrizi/repos/cloudtracker')
import cloudtracker
from cloudtracker.utility_functions import index_to_zyx, zyx_to_index
import cPickle
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
    #a list of dictionaries; each dictionary corresponds to
    #the set of all cloud shells at a given time;
    #each value in the dictionary corresponds to a list of
    #shells for a given cloud (innermost shells listed first)
    cloud_shells_list = []
    #number of shells to keep track of for each cloud
    d = 3

    LCL = 500.
    LCL_index = int(nz*dz/LCL)

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
    for clusters in clusters_list:
        #iterate over each cluster at the current time
        #increment the lifetime if it has a condensed region,
        # and find the area of xy cloud projection
        for id, cluster in clusters.iteritems():
            if cluster.has_condensed():
                lifetimes[id] = lifetimes[id] + dt
            for index in cluster.condensed_mask():
                xy_proj = set()
                proj = index_to_xy(index, ny, nx)
                xy_proj.add(proj)
                areas[id] = areas[id] + len(xy_proj)

    #filter out noise (i.e. clusters that exist only for a single time step)
    cloud_ids = np.where(lifetimes > dt)[0]

    #find cloud shells at each time step (only for clouds that exist for more than one time step)
    for clusters in clusters_list:
        #keep track of cloud shells for each cloud at the current time
        cloud_shells = {}
        for id, cluster in clusters.iteritems():
            if id in cloud_ids:
                #initialize list of shells for a given cloud
                cloud_shells[id] = []
                indices = cluster.condensed_mask()
                for x in range(d):
                    shell, expanded_indices = find_shell(indices, MC)
                    cloud_shells[id].append(shell)
                    indices = expanded_indices
        cloud_shells_list.append(cloud_shells)

   
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

    t = 20

    print "clouds at time step %2d: ", cloud_shells_list[t].keys()


    #the following is to check whether find_shell is working
    #as expected
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')

    cloud_id_small = 1645
    cloud_id_large = 2412
    
    cloud_shells = cloud_shells_list[t][cloud_id_large]

    cloud = clusters_list[t][cloud_id_large].condensed_mask()
    
    shell_one = cloud_shells[0]
    shell_points = index_to_zyx(shell_one, MC)
    x = shell_points[2,:]
    y = shell_points[1,:]
    z = shell_points[0,:]
    ax.scatter(x, y, z, c = 'y')
    shell_two = cloud_shells[1]
    shell_points = index_to_zyx(shell_two, MC)
    x = shell_points[2,:]
    y = shell_points[1,:]
    z = shell_points[0,:]
    ax.scatter(x, y, z, c = 'b')
    cloud_points = index_to_zyx(cloud, MC)
    ax.scatter(cloud_points[2,:], cloud_points[1,:], cloud_points[0,:], c = 'r')
    
    
    
    bins = np.arange(31)
    plt.figure(2)
    plt.hist(lifetimes, bins)
    plt.xlabel('lifetime (min)')
    plt.ylabel('number of clouds')
    plt.figure(3)
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

finds the shell surrounding a region in space

a shell is the set of points that are horizontally adjacent and outside the region

based off of find_halo in cloudtracker.utility_functions

inputs: indices - indices specifying the region in space
        MC - dictionary containing grid parameters (specified in model_config.cfg)

output: a numpy array with the shell indices
        a numpy array with both the closed region and shell indices

"""

def find_shell(indices, MC):

    nx = MC['nx']
    ny = MC['ny']

    #convert from indices to points (z, y, x)
    #points are ordered columnwise 
    points = index_to_zyx(indices, MC)
    stack_list = [points, ]
    #find horizontally adjacent points
    for item in ((0, -1, 0), (0, 1, 0),
                 (0, 0, -1), (0, 0, 1)):
        stack_list.append(points + np.array(item)[:, np.newaxis])

    expanded_indices = np.hstack(stack_list)

    #make sure points dont expand out of the domain
    #(wrap them if they do)
    expanded_indices[1, :] = expanded_indices[1, :]%ny
    expanded_indices[2, :] = expanded_indices[2, :]%nx

    #convert back to indices
    expanded_indices = zyx_to_index(expanded_indices[0, :],
                                    expanded_indices[1, :],
                                    expanded_indices[2, :],
                                    MC)

    expanded_indices = np.unique(expanded_indices)
    
    #only keep the indices that not in the specified region
    shell = np.setdiff1d(expanded_indices, indices, assume_unique=True)

    return shell, expanded_indices

"""
function that finds the cloud boundary at height index h

inputs: cluster - Cluster object
             h  - height index (assuming height indices start at zero)

output: array of indices corresponding to the condensed region boundary at height index h

"""

def find_boundary(cluster, h):
    nx = cluster.MC['nx']
    ny = cluster.MC['ny']
    
    boundary_3D = cluster.condensed_halo()

    hit = np.logical_and(boundary_3D >= h*nx*ny, boundary_3D <= (h+1)*nx*ny)
    return boundary_3D[hit]


if __name__ == '__main__':
    main()






    
    
    







