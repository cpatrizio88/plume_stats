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

tr = tracer_var[:]
w = w_var[:]
#qn = qn_var[:]

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






 




 
 
