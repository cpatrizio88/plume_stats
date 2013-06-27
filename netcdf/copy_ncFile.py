"""
   function that clones a ncfile (or user specified variables only)
   input the file name (file must be in working directory) and varnames

   varnames is a list of user specified variables to copy
   if varnames is empty or None then all variables will be copied

   outputs a copy of the file

"""

from netCDF4 import Dataset
import os


def copy_ncFile(filename_in, varnames):

    filename_out = filename_in + '_sample.nc'
    infile  = Dataset(filename_in+'.nc', 'r')
    #an error is thrown if the output file already exists, delete it first if it exists
    try:
        filepath = os.path.realpath(filename_out)
        os.remove(filepath)
        outfile  = Dataset(filename_out, mode='w', format='NETCDF4')
    except OSError:
        outfile  = Dataset(filename_out, mode='w', format='NETCDF4')

    print 'file', filename_out, 'created'

    #first transfer the dimensions
    for (name,value) in infile.dimensions.items():
        outfile.createDimension(name,len(value))

    #next tranfer the global attributes
    for globalatt in infile.ncattrs():
        print "copying global attribute: ",globalatt
        outfile.setncattr(globalatt,infile.getncattr(globalatt))

    #if specified by user, transfer the variable and its attributes
    #(transfer independent variables regardless)
    #otherwise transfer all variables
    #...sloppy, probably better way to do this
    if varnames:
        for (inname,invalue) in infile.variables.items():
            if inname in varnames or inname in ['x', 'y', 'z', 'time']:
                print "copying variable: ",inname
                #create the variable
                outvar = outfile.createVariable(inname,invalue.dtype,invalue.dimensions)
                #transfer the attributes
                for varattr in invalue.ncattrs():
                    outvar.setncattr(varattr,invalue.getncattr(varattr))
                #finally, copy the data
                outvar[...]=invalue[...]
    else:
         for (inname,invalue) in infile.variables.items():
            print "copying variable: ", inname
            outvar = outfile.createVariable(inname,invalue.dtype,invalue.dimensions)
            for varattr in invalue.ncattrs():
                outvar.setncattr(varattr,invalue.getncattr(varattr))
            outvar[...]=invalue[...]
         
    infile.close()

    return outfile
