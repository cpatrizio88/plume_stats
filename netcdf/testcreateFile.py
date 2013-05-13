from netCDF4 import Dataset
import numpy as np
from initialize_ncFile import initialize_ncFile
import glob

#script to combine a series of time steps of a variable of interest
#(contained in separate netCDF files) into a single netCDF file 

#-assuming the variable of interest has 3D spatial dependence


#get a list of all the netcdf BOMEX files in the current directory
ncFile_names = glob.glob('BOMEX*nc')
print ncFile_names
#use the first file to set up the dimensions of the output file
nc_in = Dataset(ncFile_names[0], 'r')
print nc_in.dimensions

xlen = len(nc_in.variables['x'])
ylen = len(nc_in.variables['y'])
zlen = len(nc_in.variables['z'])
#assume each BOMEX file corresponds to a single time step
tlen = len(ncFile_names)

#only keep track of the water vapour variable
varname = 'QV'
initialize_ncFile('testBOMEXfile.nc', varname, xlen, ylen, zlen, tlen)

nc_out = Dataset('testncfile.nc', mode = 'r+')

nc_out.variables['x'][:] = nc_in.variables['x'][:]
nc_out.variables['y'][:] = nc_in.variables['y'][:]
nc_out.variables['z'][:] = nc_in.variables['z'][:]

time = nc_out.variables['time']
#assume some value for the time step
dt = 1
tstart = 0
tend = tstart + dt*tlen
time[:] = np.arange(tstart, tend, tlen)

#loop through all of the files
#at each iteration, fill in the value of the variable of interest at the current 
#time step
for i, fname in enumerate(ncFile_names):
	nc_in = Dataset(fname, 'r')
	var = nc_in.variables[varname][...]
    #eliminate dimensions with length 1
	var = var.squeeze()
	nc_out[i, :, : ,:] = var




















