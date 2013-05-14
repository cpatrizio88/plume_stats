from netCDF4 import Dataset
import numpy as np
from initialize_ncFile import initialize_ncFile
import glob

#script to combine a series of time steps of a variable of interest
#(contained in separate netCDF files) into a single netCDF file 

#assuming the variable of interest has 3D spatial dependence

#get a list of all the netcdf BOMEX files in the current directory

ncFile_names = glob.glob('BOMEX*nc')
print "files to stack: "
for f in ncFile_names:
	print f
#use the first file to set up the dimensions of the output file
nc_in = Dataset(ncFile_names[0], 'r')

xlen = len(nc_in.variables['x'])
ylen = len(nc_in.variables['y'])
zlen = len(nc_in.variables['z'])
#assume each BOMEX file corresponds to a single time step
tlen = len(ncFile_names)

print "x length", xlen
print "y length", ylen
print "z length", zlen
print "t length", tlen

#variable to stack
varname = 'QN'
outfile = 'testBOMEXfile.nc'

print 'creating new file: ', outfile
nc_out = initialize_ncFile(outfile, varname, xlen, ylen, zlen, tlen)
nc_out.variables['x'][:] = nc_in.variables['x'][:]
nc_out.variables['y'][:] = nc_in.variables['y'][:]
nc_out.variables['z'][:] = nc_in.variables['z'][:]

time = nc_out.variables['time']
#assume some value for the time step
dt = 1
tstart = 0
tend = tstart + dt*tlen
time[:] = np.arange(tstart, tend, dt)

print 'file created, start stacking'
#loop through all of the files
#at each iteration, fill in the value of the variable of interest at the current 
#time step
for i, fname in enumerate(ncFile_names):
	print '		loading file: ', fname
	nc_in = Dataset(fname, 'r')
	varin = nc_in.variables[varname][...]
    #eliminate dimensions with length 1
	varin = varin.squeeze()
	print '		copying variable...'
	varout = nc_out.variables[varname]
	varout[i, :, : ,:] = varin

print "done stacking, the variables are: "
for var in nc_out.variables.values():
	print var




















