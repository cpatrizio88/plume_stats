from netCDF4 import Dataset
import numpy as np

nc_out = Dataset('testfile.nc', mode='w', format='NETCDF3_CLASSIC')
nc_out.createDimension('x', 10)
nc_out.createDimension('y', 10)
nc_out.createDimension('z', 10)
nc_out.createDimension('time', None)
    
x = nc_out.createVariable('x', 'f4', ('x',))
y = nc_out.createVariable('y', 'f4', ('y',))
z = nc_out.createVariable('z', 'f4', ('z',))
time = nc_out.createVariable('time', 'f4', ('time',))
    
val = nc_out.createVariable('value', 'f4', ('time', 'z', 'y', 'x'))
    
x.units = 'meters'
y.units = 'meters'
z.units = 'meters'
time.units = 'seconds since 2013-05-08 00:00:00 +0:00'
    

x[:] = np.linspace(0,9, 10)
y[:] = np.linspace(0,9, 10)
z[:] = np.linspace(0,9, 10)
time[:] = np.array([0])
val[...] = np.ones((len(time), len(z), len(y), len(x)))

for var in nc_out.variables.values():
    print var

nc_out.close()

