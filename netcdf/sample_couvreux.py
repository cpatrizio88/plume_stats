from netCDF4 import Dataset
import numpy as np
from copy_ncFile import copy_ncFile

"""
   function that samples a tracer field in a specified netcdf file 
   (there must be a vertical velocity field present in the input netcdf file as well, liquid water is optional)
   according to sampling criteria given in Crouveux (2010)

   note: the fields should be functions of (t, z, y, x)

   inputs: filename - name of netcdf file to sample from
           zb - cloud base, None if no clouds
           zt - cloud top, None if no clouds

   output:

          a numpy array of booleans, True for grid points that pass the sampling criteria, False otherwise
          (the array with be the same dimension as the fields in the netcdf file)

"""
def sample_couvreux(filename, zb, zt):

        #open the input file
    nc_in = Dataset(filename, 'r')

    tracer_var = nc_in.variables['TR01']
    w_var = nc_in.variables['W']
    z_var = nc_in.variables['z']

    tr = tracer_var[:]
    w = w_var[:]
    z = z_var[:]

    #if there are clouds, load liquid water variable
    if zb:
       qn_var = nc_in.variables['QN']
       qn = qn_var[:]

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
    #(according to Couvreux, 2009) at a given time and height index
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

    #find points that pass sampling criteria
    print 'sampling...'
    hit = np.logical_and(tr > tr_threshold, w > 0)
    #if there are clouds, add condition liquid water content > 0
    #to sampling criteria above height z - zb + (zt - zb)/4
    if zb:
            #find height index corresponding to z =  zb + (zt - zb)/4
            cld_index = np.where(z > zb + (zt - zb)/4)[0][0]
            #update points that will be sampled
            hit[:,cld_index::,:,:] = np.logical_and(hit[:,cld_index::,:,:], qn[:,cld_index::,:,:] > 0)

    #remove single dimensions from hit
    #(e.g. when sampling a single time slice, the netcdf file typically specifies time to be dimension of length 1)
    hit = np.squeeze(hit)

    nc_in.close()

    return hit

        






















 




 
 
