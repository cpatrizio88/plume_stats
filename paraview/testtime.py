from netCDF4 import Dataset
import site
site.addsitedir('C:\\Users\\Casey\\atscCode\\A405\\python\\thermlib')
site.addsitedir('C:\\Users\\Casey\\atscCode\\A405\\python\\skew_T')
site.addsitedir('C:\\Users\\Casey\\atscCode\\A405\\python\\ode45')
from scipy.integrate import ode
import numpy as np
from constants import constants as c
from new_thermo import thetaep, tinvert_thetae, wsat
from calcBuoy import calcBuoy
from nudge import nudge

def testtime():
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
    wTcloud = wsat(dewpoint[p900_level] + c.Tc, press[p900_level]*100.)
    thetaeVal=thetaep(dewpoint[p900_level] + c.Tc, temp[p900_level] + c.Tc, press[p900_level]*100.)
    
    height_800=height[p800_level]
    
    Yinit = [0.5, height_800]  #(intial velocity = 0.5 m/s, initial height in m)
    tinit = 0
    tfin = 1000
    dt = 5
    
    #want to integrate F using ode45 (from MATLAB) equivalent integrator
    r = ode(F).set_integrator('dopri5')
    r.set_f_params(thetaeVal, interpTenv, interpTdEnv, interpPress)
    r.set_initial_value(Yinit, tinit)
    
    Y = np.array(Yinit)
    t = np.array(tinit)
    
    #stop integration when the parcel changes direction, or time runs out
    #while r.successful() and r.t < tfin and r.y[0] > 0:
    while r.successful() and r.t < tfin: 
        #find y at the next time step
        #(r.integrate(t) updates the fields r.y and r.t so that r.Y = F(t) and r.t = t 
        #where F is the function being integrated)
        r.integrate(r.t+dt)
        #keep track of Y (velocity and position) at each time step
        Y = np.vstack((Y, r.y))
        t = np.vstack((t, r.t))

    #create netCDF file to hold time evolution of cloud parcel
    nc_out = Dataset('adiabatic_ascent.nc', 'w', format = 'NETCDF3_CLASSIC')
    nc_out.createDimension('x', 10)
    nc_out.createDimension('y', 10)
    nc_out.createDimension('z', 200)
    nc_out.createDimension('time', len(t))
    
    x = nc_out.createVariable('x', 'f4', ('x',))
    y = nc_out.createVariable('y', 'f4', ('y',))
    z = nc_out.createVariable('z', 'f4', ('z',))
    time = nc_out.createVariable('time', 'f4', ('time',))
    
    w_l = nc_out.createVariable('liquid_water', 'f4', ('time', 'z', 'y', 'x'))
    
    x.units = 'meters'
    y.units = 'meters'
    z.units = 'meters'
    time.units = 'seconds since 2013-05-08 00:00:00 +0:00'
    w_l.units = 'kg liquid water per kg air'
    
    x[:] = np.linspace(-5e3, 5e3, len(x))
    y[:] = np.linspace(-5e3, 5e3, len(y))
    z[:] = np.linspace(height[0], height[-1], len(z))
    time[:] = t   
    
    cloud_height = Y[:,1]   
    Tcloud = np.zeros(cloud_height.size)
    wvCloud = np.zeros(cloud_height.size)
    wlCloud = np.zeros(cloud_height.size)
    
    for i in range(0, len(t)):
        the_press = interpPress(cloud_height[i])*100.
        Tcloud[i], wvCloud[i], wlCloud[i] = tinvert_thetae(thetaeVal,
                                                        wTcloud, the_press)
        print wlCloud[i]
                                                           
       
    #fill in the liquid water field at each time step  
    for tindex, tval in enumerate(t):
        #first initialize to zero
        w_l[tindex,:,:,:] = np.zeros((len(x), len(y), len(z)))
        #find index in z closest to the current cloud height
        cloud_height_index = np.where(abs(cloud_height[tindex] - z) < 100.)[0][0]
        #fill in the liquid water content at this height
        w_l[tindex,cloud_height_index,:,:] = wlCloud[tindex]*np.ones((len(y), len(x)))
            
#F returns the buoyancy (and height) at a given time step and height
def F(t, y, thetae0, interpTenv, interpTdEnv, interpPress):
    #y[0] is the velocity, y[1] is the height
    yp = np.zeros((2,1))
    #yp[0] is the buoyancy (acceleration), yp[1] is the velocity
    yp[0] = calcBuoy(y[1], thetae0, interpTenv, interpTdEnv, interpPress)
    yp[1] = y[0] 
    return yp
    
    
if __name__ == "__main__":
    testtime()
    
    
    
    