from netCDF4 import Dataset
import numpy as np
from copy_ncFile import copy_ncFile

#name of nc file to sample from (must be in working directory)
filename_tr = 'testBOMEXtracer_expdecay'

#set up output nc file
nc_out = copy_ncFile(filename_tr)

#nc_w = Dataset('testBOMEXtracer_w.nc', 'r')
#nc_qn = Dataset('testBOMEXtracer_QN.nc',  'r')

tracer_var = nc_out.variables['TR01']
#w_var = nc_w.variables['W']
#qn_var = nc_qn.variables['QN']
z_var = nc_out.variables['z']


tr = tracer_var[:]
#w = w_var[:]
#qn = qn_var[:]
z = z_var[:]

dz = np.diff(z)
#dz is has one less element than z, append an element (average height spacing) to the end
dz = np.append(dz, np.mean(dz))

tlen = tr.shape[0]
print "tlen", tlen
zlen = tr.shape[1]
print "zlen", zlen
xlen = tr.shape[2]
ylen = tr.shape[3]
#compute standard deviation and mean at each height level, first flatten x-y domain
#tr_stdev and tr_mean should have shape (tlen, zlen)
tr_flat = tr.reshape(tlen, zlen, xlen*ylen)
tr_stdev = np.std(tr_flat, axis=2)
print "tr_stdev.shape", tr_stdev.shape
tr_mean = np.mean(tr_flat, axis=2)
print "tr_mean.shape", tr_mean.shape

#function to compute the minimum standard deviation 
#(according to Couvreux, 2009) at a given time and height index
stdev_min = lambda i, j: (0.05/z[j])*sum(tr_stdev[i, 0:j]*dz[0:j])

tr_sample = np.zeros(tr.shape)
#scaling factor for sampling threshold
m = 1
print 'sampling...'
for i in np.arange(tlen):
	print "time index: ", i
	for j in np.arange(zlen):
		print "height index: ", j
		tr_threshold = tr_mean[i, j] + m*max(tr_stdev[i, j], stdev_min(i, j)) 
	    #initialize logical indexing array to be the same shape as tracer array 
		#(w/ all False entries)
		hit = np.zeros(tr.shape) > 1
		#find the grid cells in the current x-y slice that should be sampled
		hit[i, j, :, :] = tr[i, j, :, :] > tr_threshold
		tr_sample[hit] = tr[hit]
#fill netcdf tracer variable with sampled points
tracer_var[:] = tr_sample
tr
nc_out.close()





















 




 
 
