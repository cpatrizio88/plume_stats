
from netCDF4 import Dataset
import numpy as np
import os

#initializes a netCDF file to hold a single variable with dimensions x, y, z, time
#dimension lengths are specified by user
def initialize_ncFile(filename, varname, xlen, ylen, zlen, tlen): 
    

    #an error is thrown if the file already exists, delete it first if it exists
    try:
        filepath = os.path.realpath(filename)
        os.remove(filepath)
        nc_out = Dataset(filename, mode='w', format='NETCDF4')
    except OSError:
        nc_out = Dataset(filename, mode='w', format='NETCDF4')
    
    nc_out.createDimension('x', xlen)
    nc_out.createDimension('y', ylen)
    nc_out.createDimension('z', zlen)
    nc_out.createDimension('time', tlen)
    
    x = nc_out.createVariable('x', 'f8', ('x',))
    y = nc_out.createVariable('y', 'f8', ('y',))
    z = nc_out.createVariable('z', 'f8', ('z',))
    time = nc_out.createVariable('time', 'f8', ('time',))
    val = nc_out.createVariable(varname, 'f8', ('time', 'z', 'y', 'x'))

    x.units = 'meters'
    y.units = 'meters'
    z.units = 'meters'
    time.units = 'seconds since 2013-05-08 00:00:00 +0:00'
   
    return nc_out
