import glob
from netCDF4 import Dataset
from sampling_funcs import sample_couvreux, sample_condensed, sample_core
from initialize_ncFile import initialize_ncFile
import numpy as np

"""
script that produces netcdf files required for cloudtracker

samples all netcdf files in the current directory and outputs
netcdf files containing core, condensed and plume variables which are arrays of 1's and 0's
(to indicate which grid points pass the given sampling criteria)

(for details see: https://github.com/freedryk/cloudtracker/blob/master/README.txt)

*note: assuming all netcdf files in the current directory correspond to single time steps

"""

nc_filenames = glob.glob('/tera/phil/cloudtracking/*[!_cldtrcksample].nc')

sample_liqwater = False

for filename in nc_filenames:
    print 'opening file: ', filename
    nc_in = Dataset(filename, 'r')
    xlen  = len(nc_in.dimensions['x'])
    ylen = len(nc_in.dimensions['y'])
    zlen = len(nc_in.dimensions['z'])
    varnames = ['core', 'condensed', 'plume', 'u', 'v', 'w']
    nc_out = initialize_ncFile(filename[:-3]+'_cldtrcksample.nc', varnames, xlen, ylen, zlen)
    #fill velocity variables
    u_var = nc_out.variables['u']
    v_var = nc_out.variables['v']
    w_var = nc_out.variables['w']
    u_var[:] = nc_in.variables['U'][:]
    v_var[:] = nc_in.variables['V'][:]
    w_var[:] = nc_in.variables['W'][:]
    nc_in.close()
    #find plume, core and condensed regions
    plume  = sample_couvreux(filename, sample_liqwater)
    #convert from True/False to 1/0
    plume = plume/1
    core = sample_core(filename, sample_liqwater)
    core = core/1
    if sample_liqwater:
        condensed = sample_condensed(filename)
        condensed = condensed/1
    else:
        condensed = np.zeros((zlen, ylen, xlen))
    #fill variable in nc_out with points to sample
    plume_var = nc_out.variables['plume']
    plume_var[:] = plume
    core_var = nc_out.variables['core']
    core_var[:] = core
    condensed_var = nc_out.variables['condensed']
    condensed_var[:] = condensed
    nc_out.close()
    
