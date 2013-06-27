from netCDF4 import Dataset
import numpy as np
from copy_ncFile import copy_ncFile

#name of nc file to sample from (must be in working directory)
#(assuming all variables are functions of t, z, y, x)
filename = 'testNCHAPP1tracer'
#open input file 
nc_in = Dataset(filename+'.nc', 'r')
#set up the output file by copying only the tracer variable from the input file
nc_out = copy_ncFile(filename, ['TR01'])

tracer_var = nc_out.variables['TR01']
w_var = nc_in.variables['W']
#qn_var = nc_qn.variables['QN']
z_var = nc_in.variables['z']


tr = tracer_var[:]
w = w_var[:]
#qn = qn_var[:]
z = z_var[:]

dz = np.diff(z)
#dz has one less element than z, append an element to the end
dz = np.append(dz, np.mean(dz))

#cloud base, None if dry case
zb = None
#cloud top
zt = None

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
tr_sample = np.zeros(tr.shape)
tr_threshold = np.zeros(tr.shape)
#scaling factor for sampling threshold
m = 1

print 'calculating tracer thresholds...'
for i in np.arange(tlen):
	for j in np.arange(zlen):
               #keep track of thresholds at each time and height
		tr_threshold[i, j, :, :] = tr_mean[i, j] + m*max(tr_stdev[i, j], stdev_min(i, j)) 
	
#fill netcdf tracer variable with sampled points
print 'sampling...'
hit = np.logical_and(tr > tr_threshold, w > 0)
tr_sample[hit] = tr[hit]
tracer_var[:] = tr_sample
#if there are clouds, add condition liquid water content > 0
#to sampling criteria above height z - zb + (zt - zb)/4
if zb:
        #find height index corresponding to z =  zb + (zt - zb)/4
        cld_index = np.where(z > zb + (zt - zb)/4)[0][0]
        #update points that will be sampled
        hit[:,cld_index::,:,:] = np.logical_and(hit[:,cld_index::,:,:], qn[:,cld_index::,:,:] > 0)
        tr_sample[hit] = tr[hit]
        tracer_var[:] = tr_sample
        
nc_out.close()
nc_in.close()






















 




 
 
