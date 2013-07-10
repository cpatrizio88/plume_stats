from netCDF4 import Dataset
import numpy as np

#Module that contains sampling functions used to produce input for cloud tracking algorithm

"""
   function that samples a tracer field in a specified netcdf file 
   (there must be a vertical velocity field present in the input netcdf file as well, liquid water is optional)
   according to sampling criteria given in Couvreux (2010)

   note: the fields should be functions of (t, z, y, x)

   inputs: filename - name of netcdf file to sample from
           optional parameter - True if adding liquid water to sampling criteria 
           

   output:

          a numpy array of booleans, True for grid points that pass the sampling criteria, False otherwise
          (the array with be the same dimension as the fields in the netcdf file)
"""
def sample_couvreux(filename, *args):
    if args:
        sample_liqwater = args[0]
    else:
        sample_liqwater = False
    #open the input file
    nc_in = Dataset(filename, 'r')

    tracer_var = nc_in.variables['TR01']
    w_var = nc_in.variables['W']
    z_var = nc_in.variables['z']

    tr = tracer_var[:]
    w = w_var[:]
    z = z_var[:]

    if sample_liqwater:
       qn_var = nc_in.variables['QN']
       qn = qn_var[:]

    tr_thresholds = compute_tracer_thresholds(tr, z)

    #find points that pass sampling criteria
    print 'sampling plume...'
    hit = np.logical_and(tr > tr_threshold, w > 0)
    if sample_liqwater:
        hit = np.logical_or(hit, qn > 0)

    #additional criteria specified in Couvreux, not using for now:
    #if there are clouds, add condition liquid water content > 0
    #to sampling criteria above height z - zb + (zt - zb)/4
    #if args:
    #        zb = args[0]
    #        zt = args[1]
    #         #find height index corresponding to z =  zb + (zt - zb)/4
    #        cld_index = np.where(z > zb + (zt - zb)/4)[0][0]
    #        #update points that will be sampled
    #        hit[:,cld_index::,:,:] = np.logical_and(hit[:,cld_index::,:,:], qn[:,cld_index::,:,:] > 0)

    #remove single dimensions from hit
    #(e.g. when sampling a single time slice from a netcdf file, the time dimension may be length 1)
    hit = np.squeeze(hit)
    nc_in.close()
    return hit

"""
  sample fields contained in a netcdf file using core sampling criteria

  input: filename - netcdf file to sample from
         optional parameter - True if adding liquid water to sampling criteria
         
  output: a numpy array of booleans, for dry case: True for positively buoyant, upward moving grid cells with tracer
                                                   concentration above a threshold, False otherwise
                                     for liquid water case: True for positively buoyant, upward moving grid cells wit                                                            with liquid water, False otherwise
                                      

"""

def sample_core(filename, *args):
    if args:
      sample_liqwater = args[0]
    else:
      sample_liqwater = False
    nc_in = Dataset(filename, 'r')
    if sample_liqwater:
        qn = nc_in.variables['QN'][:]
        qn = np.squeeze(qn)
    T = nc_in.variables['TABS'][:]
    T = np.squeeze(T)
    qv = nc_in.variables['QV'][:]
    qv = np.squeeze(qv)
    p = nc_in.variables['P'][:]
    p = np.squeeze(p)
    w = nc_in.variables['W'][:]
    w = np.squeeze(w)
    tr = nc_in.variables['TR01'][:]
    z = nc_in.variables['z'][:]

    #calculate thetav from T, p, and qv
    p0 = 1000*100 #reference pressure (MSL), Pa
    Rd = 287. #gas constant for dry air, J/kg/K
    Rv = 461. #gas constant for water vapour, J/kg/K
    cp = 1004 #specific heat of air J/kg/K
    eps = Rd/Rv
    Tv = T*(1+qv/eps)/(1+qv)
    thetav = Tv*(p0/p)**(Rd/cp)
    #calculate delta_thetav, the thetav anomaly (difference between thetav and mean thetav at a given height)
    #this is proportional to the buoyancy
    print 'sampling core...'
    thetav_z = np.mean(np.mean(thetav, axis=1), axis=1)
    zlen = len(thetav_z)
    thetav_z = thetav_z.reshape((zlen,1,1))
    delta_thetav = thetav - thetav_z
    #find positively buoyant, upward moving grid cells 
    hit = np.logical_and(delta_thetav > 0, w > 0)
    if sample_liqwater:
        hit = np.logical_and(qn > 0, hit)
        return hit
    else:
        tr_thresholds = compute_tracer_thresholds(tr, z)
        hit = np.logical_and(tr > tr_thresholds, hit)
        return hit
        


def sample_condensed(filename):
    print 'sampling condensed...'
    nc_in = Dataset(filename)
    qn = nc_in.variables['QN']
    hit = qn > 0
    hit = np.squeeze(hit)
    return hit

"""
returns tracer thresholds (according to Couvreux, 2010) at height levels specified in z

   inputs: tr - tracer concentration field (must be a function of t, z, y, x)
           z  - heights (m)

   output: array containing tracer thresholds (same dimension as tr)

"""

def compute_tracer_thresholds(tr, z):

    dz = np.diff(z)
    #dz has one less element than z, append an element to the end
    dz = np.append(dz, np.mean(dz))

    tlen = tr.shape[0]
    #print "tlen", tlen
    zlen = tr.shape[1]
    #print "zlen", zlen
    xlen = tr.shape[2]
    ylen = tr.shape[3]
    #compute standard deviation and mean at each height level, first flatten x-y domain
    #tr_stdev and tr_mean should have shape (tlen, zlen)
    tr_flat = tr.reshape(tlen, zlen, xlen*ylen)
    tr_stdev = np.std(tr_flat, axis=2)
    #print "tr_stdev.shape", tr_stdev.shape
    tr_mean = np.mean(tr_flat, axis=2)
    #print "tr_mean.shape", tr_mean.shape

    #function to compute the minimum standard deviation 
    #(according to Couvreux, 2010) at a given time and height index
    stdev_min = lambda i, j: (0.05/z[j])*sum(tr_stdev[i, 0:j]*dz[0:j])

    #initialize arrays
    tr_threshold = np.zeros(tr.shape)
    #scaling factor for sampling threshold
    m = 1
    print 'calculating tracer thresholds...'
    for i in np.arange(tlen):
            for j in np.arange(zlen):
                    #keep track of thresholds at each time and height
                            tr_threshold[i, j, :, :] = tr_mean[i, j] + m*max(tr_stdev[i, j], stdev_min(i, j))

    return tr_threshold
    
