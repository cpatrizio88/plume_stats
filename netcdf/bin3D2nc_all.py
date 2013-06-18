import os
import glob

# converts all bin3D files in the current directory to netcdf files 

binFile_names = glob.glob('*.bin3D')
#path to bin3D2nc script, which converts a single bin3D file to .nc file
fn = '/tera/cpatrizi/test2/bin/bin3D2nc'


for f in binFile_names:
    cmd = fn + ' ' + f
    os.system(cmd)
    
