from netCDF4 import Dataset
import numpy as np
from initialize_ncFile import initialize_ncFile
import glob

#script to combine a series of time steps of variables
#(the time steps are contained in separate netCDF files) into a single netCDF file 

#assuming all variables have 3D spatial dependence

#get a list of all the netcdf files in the current directory
ncFile_names = glob.glob('NCHAPP1*nc')
#sort the files chronologically (assuming they are time stamped)
ncFile_names.sort()
print "files to stack: "
for f in ncFile_names:
	print f
#use the first file to set up the dimensions of the output file
nc_in = Dataset(ncFile_names[0], 'r')
xlen = len(nc_in.variables['x'])
ylen = len(nc_in.variables['y'])
zlen = len(nc_in.variables['z'])
#assume each netcdf file corresponds to a single time step
tlen = len(ncFile_names)
print "x length", xlen
print "y length", ylen
print "z length", zlen
print "t length", tlen
#variables to stack
varnames = ['TR01', 'W']
outfile = 'testNCHAPP1tracer.nc'
print 'creating new file: ', outfile
nc_out = initialize_ncFile(outfile, varnames,  xlen, ylen, zlen, tlen)
#set up units, name, etc
for varname in varnames:
        varin = nc_in.variables[varname]
        varout = nc_out.variables[varname]
        for varattr in varin.ncattrs():
                varout.setncattr(varattr,varin.getncattr(varattr))
#fill in dimensions
x = nc_out.variables['x']
y = nc_out.variables['y']
z = nc_out.variables['z']
x[:] = nc_in.variables['x'][:]
y[:] = nc_in.variables['y'][:]
z[:] = nc_in.variables['z'][:]
time = nc_out.variables['time']
#the time step is 4 min. for this particular case
dt = 2
tstart = 0
tend = tstart + dt*tlen
time[:] = np.arange(tstart, tend, dt)
time.units = 'minutes since 2013-05-08 00:00:00 +0:00'
print 'file created, start stacking'

#loop through all of the files, at each iteration
#copy over the values of the variables
print 'copying variables' , varnames
for i, fname in enumerate(ncFile_names):
	print '		loading file: ', fname
	nc_in = Dataset(fname, 'r')
        for varname in varnames:
            varin = nc_in.variables[varname][...]
            #eliminate dimensions with length 1
            varin = varin.squeeze()
            print '		copying variable', varname
            varout = nc_out.variables[varname]
            varout[i, :, : ,:] = varin

print "done stacking, the variables are: "
for var in nc_out.variables.values():
	print var

nc_out.close()
nc_in.close()




















