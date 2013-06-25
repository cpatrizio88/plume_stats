from netCDF4 import Dataset
import numpy as np
from copy_ncFile import copy_ncFile

#name of nc file to sample from (must be in working directory)
filename_tr = 'testBOMEXtracer_expdecay'

#set up output nc file
nc_out = copy_ncFile(filename_tr)

nc_w = Dataset('testBOMEXtracer_w.nc', 'r')
#nc_qn = Dataset('testBOMEXtracer_QN.nc',  'r')

tracer_var = nc_out.variables['TR01']
w_var = nc_w.variables['W']
#qn_var = nc_qn.variables['QN']
z_var = nc_out.variables['z']


tr = tracer_var[:]
w = w_var[:]
#qn = qn_var[:]
z = z_var[:]
dz = np.diff(z)

#not sure about this line, want to compute standard deviation
#at each height level
#first need to flatten the xy domain 
tlen = tr.shape[0]
zlen = tr.shape[1]
xlen = tr.shape[2]
ylen = tr.shape[3]
tr_stdev = np.std(tr.reshape(tlen, zlen, xlen*ylen), axis=2)
#tr_stdev should have shape (tlen, zlen)

#function to compute the minimum standard deviation 
#(according to Couvreux, 2009) at a given time and height
stdev_minFn = lambda i, h :: (0.05/h)*sum(tr_stdev[i, 0:z.where(h)]*dz[0:z.where(h)])

stdev_min = stdev_minFn(arange(0, tlen), z)

m = 1

tracer_new = np.zeros(tr.shape)

for i in arange(tlen):
	for j in arange(zlen):
		tr_threshold = m*max(tr_stdev(i, j), stdev_min(i, j))
		hit = tr[i, j, :, :] > tr_threshold
		tracer_new[hit] = tr[hit]













tr_threshold = 1
w_threshold = 0
#qn_threshold = 0

#logical_and compares truth values element-wise
print 'sampling...'
hit = np.logical_and(tr > tr_threshold, w > w_threshold)
#hit  = np.logical_and(np.logical_and(tr > tr_threshold,  w > w_threshold), qn > qn_threshold)
tracer_new = np.zeros(tr.shape)
tracer_new[hit] = tr[hit]
#update tracer variable with sampled points
tracer_var[:] = tracer_new

nc_out.close()
nc_w.close()









 




 
 
