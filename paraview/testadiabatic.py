from netCDF4 import Dataset
import site
site.addsitedir('C:\\Users\\Casey\\atsc code\\A405\\python\\thermlib')
from scipy.integrate import ode
import numpy as np
from constants import constants as c
from nudge import nudge

def testadiabatic():
    #create netCDF file to hold data for adiabatic ascent of air parcel 
    #(using littlerock sounding)
    
    nc_out = Dataset('adiabatic_ascent.nc', 'w')
   
    nc_out.createDimension('z', None)
    nc_out.createDimension('time', None)
    
    
    
    
    
    filename = 'littlerock.nc'
    print 'reading file: %s\n' %(filename)
    nc_file = Dataset(filename)
    var_names = nc_file.variables.keys()
    print nc_file.ncattrs()
    print nc_file.units
    print nc_file.col_names
    
    sound_var = nc_file.variables[var_names[3]]
    press = sound_var[:,0]
    height = sound_var[:,1]
    temp = sound_var[:,2]
    dewpoint = sound_var[:,3]
    
    #height must have unique values
    newHeight= nudge(height)
    #Tenv and TdEnv interpolators return temp. in deg C, given height in m
    #Press interpolator returns pressure in hPa given height in m
    interpTenv = lambda zVals: np.interp(zVals, newHeight, temp)
    interpTdEnv = lambda zVals: np.interp(zVals, newHeight, dewpoint)
    interpPress = lambda zVals: np.interp(zVals, newHeight, press)
    p900_level = np.where(abs(900 - press) < 2.)
    p800_level = np.where(abs(800 - press) < 7.)
    thetaeVal=thetaep(dewpoint[p900_level] + c.Tc,temp[p900_level] + c.Tc,press[p900_level]*100.)
    
    height_800=height[p800_level]
    
    yinit = [0.5, height_800]  #(intial velocity = 0.5 m/s, initial height in m)
    tinit = 0
    tfin = 2500
    dt = 10
    
    #want to integrate F using ode45 (from MATLAB) equivalent integrator
    r = ode(F).set_integrator('dopri5')
    r.set_f_params(thetaeVal, interpTenv, interpTdEnv, interpPress)
    r.set_initial_value(yinit, tinit)
    
    y = np.array(yinit)
    t = np.array(tinit)
    
    #stop integration when the parcel changes direction, or time runs out
    while r.successful() and r.t < tfin and r.y[0] > 0:
        #find y at the next time step
        #(r.integrate(t) updates the fields r.y and r.t so that r.y = F(t) and r.t = t 
        #where F is the function being integrated)
        r.integrate(r.t+dt)
        #keep track of y at each time step
        y = np.vstack((y, r.y))
        t = np.vstack((t, r.t))
        
cloud_height=y[:,1]
print len(cloud_height)