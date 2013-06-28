
from netCDF4 import Dataset
import numpy as np
import os

"""
initializes a netCDF file to hold multiple variables with dimensions x, y, z, time

inputs: xlen, ylen, zlen - spatial dimension lengths
        varnames - list with variable names
        optional argument - specify time dimension length if needed

"""
def initialize_ncFile(filename, varnames, xlen, ylen, zlen, *args):
    

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
    x = nc_out.createVariable('x', 'f8', ('x',))
    y = nc_out.createVariable('y', 'f8', ('y',))
    z = nc_out.createVariable('z', 'f8', ('z',))
    
    if args:
        tlen = args[0]
        nc_out.createDimension('time', tlen)
        time = nc_out.createVariable('time', 'f8', ('time',))
        for varname in varnames:
            val =  nc_out.createVariable(varname, 'f8', ('time', 'z', 'y', 'x'))
    else:
        for varname in varnames:
            val =  nc_out.createVariable(varname, 'f8', ('z', 'y', 'x'))
        
    x.units = 'meters'
    y.units = 'meters'
    z.units = 'meters'
    if tlen:
        time.units = 'seconds since 2013-05-08 00:00:00 +0:00'
   
    return nc_out
