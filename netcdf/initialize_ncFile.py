from netCDF4 import Dataset
import numpy as np

#creates a netCDF file to hold a single variable with dimensions x, y, z, time
def initialize_ncFile(filename, varname, xlen, ylen, zlen, tlen): 

	nc_out = Dataset(filename, mode='w', format='NETCDF3_CLASSIC')
	nc_out.createDimension('x', xlen)
	nc_out.createDimension('y', ylen)
	nc_out.createDimension('z', zlen)
	nc_out.createDimension('time', tlen)

	x = nc_out.createVariable('x', 'f4', ('x',))
	y = nc_out.createVariable('y', 'f4', ('y',))
	z = nc_out.createVariable('z', 'f4', ('z',))
	time = nc_out.createVariable('time', 'f4', ('time',))

	val = nc_out.createVariable(varname, 'f4', ('time', 'z', 'y', 'x'))

	x.units = 'meters'
	y.units = 'meters'
	z.units = 'meters'
	time.units = 'seconds since 2013-05-08 00:00:00 +0:00'

	nc_out.close()
