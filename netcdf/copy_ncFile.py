"""function that clones a ncfile
   input the file name (file must be in working directory)
   outputs a copy of the file
"""

from netCDF4 import Dataset
import os

def copy_ncFile(filename_in):

    filename_out = filename_in + '_sample.nc'
    infile  = Dataset(filename_in+'.nc', 'r')
    #an error is thrown if the output file already exists, delete it first if it exists
    try:
        filepath = os.path.realpath(filename_out)
        os.remove(filepath)
        outfile  = Dataset(filename_out, mode='w', format='NETCDF4')
    except OSError:
        outfile  = Dataset(filename_out, mode='w', format='NETCDF4')
    

    #first transfer the dimensions
    for (name,value) in infile.dimensions.items():
        outfile.createDimension(name,len(value))

    #next tranfer the global attributes
    for globalatt in infile.ncattrs():
        print "copying global attribute: ",globalatt
        outfile.setncattr(globalatt,infile.getncattr(globalatt))

    #transfer each variable and its attributes
    for (inname,invalue) in infile.variables.items():
        print "copying variable: ",inname
        #create the variable
        outvar = outfile.createVariable(inname,invalue.dtype,invalue.dimensions)
        #transfer the attributes
        for varattr in invalue.ncattrs():
            outvar.setncattr(varattr,invalue.getncattr(varattr))
        ## #finally, copy the data
        outvar[...]=invalue[...]
         
        
    infile.close()

    
    return outfile
