from netCDF4 import Dataset
import numpy as np
#import os



soundings = Dataset('littlerock.nc', 'r')

print soundings.variables

print soundings.col_names

var_names = soundings.variables.keys()

sound_var = soundings.variables[var_names[3]]
height = sound_var[:,1]
temp_data = sound_var[:,2]

#nc_out = Dataset('testsounding.nc', 'w', clobber = True, format='NETCDF4')
#os.remove('C:\\Users\\Casey\\repos\\casey\\paraview\\testsounding.nc')

nc_out = Dataset('testsounding.nc', 'w', clobber = True, format = 'NETCDF4')

nc_out.createDimension('x', 50)
nc_out.createDimension('y', 50)
nc_out.createDimension('z', len(height))

x = nc_out.createVariable('x', 'f4', ('x',))
y = nc_out.createVariable('y', 'f4', ('y',))
z = nc_out.createVariable('z', 'f4', ('z',))

x.units = 'meters'
y.units = 'meters'
z.units = 'meters'

x[:] = np.linspace(-1000, 1000, len(x))
y[:] = np.linspace(-1000, 1000, len(y))
z[:] = np.linspace(height[0], height[-1], len(temp_data))

temp_field = nc_out.createVariable('temp_field', 'f4', ('x', 'y', 'z'))

temp_field.units = 'deg C'

print nc_out.variables.keys()
temp_interp = np.interp(z, height, temp_data)


for k, zval in enumerate(z):
    temp_field[:,:,k] = temp_interp[k]*np.ones((len(x), len(y)))


nc_out.close()
soundings.close()
