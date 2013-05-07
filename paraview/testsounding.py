from netCDF4 import Dataset
import numpy as np

soundings = Dataset('littlerock.nc', 'r')

print soundings.variables

print soundings.col_names

var_names = soundings.variables.keys()

sound_var = soundings.variables[var_names[3]]


nc_out = Dataset('testsounding.nc', 'w', clobber = True, format='NETCDF4')



nc_out.createDimension('x', 50)
nc_out.createDimension('y', 50)
nc_out.createDimension('z', None)
#nc_out.createDimension('time', None)

x = nc_out.createVariable('x', 'f4', ('x',))
y = nc_out.createVariable('y', 'f4', ('y',))
z = nc_out.createVariable('z', 'f4', ('z',))
#time = nc_out.createVariable('time', 'f4', ('time',))

x = np.arange(-1000, 1000, 100)
y = np.arange(-1000, 1000, 100)

height = sound_var[:,1]
z = height

temp = sound_var[:,2]

temp_field = nc_out.createVariable('temp_field', 'f4', ('x', 'y', 'z'))

for i in x:
    for j in y:
        temp[i,j,:] = temp




nc_out.close()
soundings.close()
