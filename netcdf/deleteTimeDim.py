from netCDF4 import Dataset
import glob

#script to remove time dimension (assuming it is length 1) 
#from variables in a netcdf file 

ncfile_in=glob.glob("[!new]*nc")[0]
print "loading file: %s" %(ncfile_in,)
ncfile_out="new{}".format(ncfile_in)
infile = Dataset(ncfile_in, 'r')
outfile = Dataset(ncfile_out, mode='w',format='NETCDF3_CLASSIC')

for (name, value) in infile.dimensions.items():
    if name != 'time':
        outfile.createDimension(name, len(value))

for globalatt in infile.ncattrs():
    print "copying global attribute: ", globalatt
    outfile.setncattr(globalatt,infile.getncattr(globalatt))
	

for (inname, invar) in infile.variables.items():
    dim = invar.dimensions
    if 'time' in dim and len(dim) != 1:
        newdim = tuple(y for y in dim if y != 'time')
        outvar = outfile.createVariable(inname, invar.dtype, newdim)
        newvar = invar[...]
        newvar = newvar.squeeze()
        outvar[:] = newvar
    elif 'time' not in dim:
        outvar = outfile.createVariable(inname, invar.dtype, dim)
        outvar[:] = invar[...]
        
    for varattr in invar.ncattrs():
        outvar.setncattr(varattr, invar.getncattr(varattr))

    for var in outfile.variables.values():
        print var


        
  
 


    

    


 
    


