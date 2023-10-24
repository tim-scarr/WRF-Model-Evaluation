from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.pyplot import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import cartopy.feature as cfeature
from mpl_toolkits import basemap
import pandas
import numpy as np
import os
import string
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, CoordPair, extract_times, ALL_TIMES)
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter


# Open the NetCDF file
ncfile = Dataset("wrfout_d01_2017-03-20_00:00:00")



# # Get the sea level pressure

slp = getvar(ncfile, "slp", timeidx=0)
lat, lon  = latlon_coords(slp)

# #Smooth sea level pressure - noisy near mountains

smooth_slp = smooth2d(slp, 3, cenweight=4)

# #Get the latitude and longitude points

lats, lons = latlon_coords(slp)

# #Cartopy mapping object
cart_proj = get_cartopy(slp)
alphabet = list(string.ascii_lowercase)
start_number = 8
array_length = 24

array_times = [start_number + 2 * i for i in range(array_length)]

time     =  extract_times(ncfile, timeidx=ALL_TIMES, method='cat', meta=False)
times    =  pandas.to_datetime(time).strftime('%Y-%m-%d_%H:%M').values



def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def wrf_vort( U, V, dx ):
    """ Calculate the relative vorticity given the U and V vector components in m/s
    and the grid spacing dx in meters. U and V must be the same shape.
    returns: numpy.ndarray of vorticity values s^-1 same shape as U and V """
    assert U.shape == V.shape, 'Arrays are different shapes. They must be the same shape.'
    dy = dx
    du = np.gradient( U )
    dv = np.gradient( V )
    return ( dv[-1]/dx - du[-2]/dy )

def cal_vor(t):
    dx        = ncfile.DX
    ## Extract ref vars
    ua        = getvar(ncfile, "ua",      timeidx=t)
    va        = getvar(ncfile, "va",      timeidx=t)
    p         = getvar(ncfile, "pressure",timeidx=t)
    u850      = to_np(interplevel(ua, p, 850))
    v850      = to_np(interplevel(va, p, 850))
    u500      = to_np(interplevel(ua, p, 500))
    v500      = to_np(interplevel(va, p, 500))
    u200      = to_np(interplevel(ua, p, 200))
    v200      = to_np(interplevel(va, p, 200))
    # Create xarrays
    rvor850   = wrf_vort(u850, v850, dx)
    rvor200   = wrf_vort(u200, v200, dx)
    rvor500   = wrf_vort(u500, v500, dx)
    lat, lon  = latlon_coords(ua)
    return rvor850,rvor200,lat, lon,rvor500


def geopotential_500(time,a):
    
    
    
    z = smooth2d(getvar(ncfile, "z", timeidx = time ), 3, cenweight=1)
    z_meta = getvar(ncfile, "z", timeidx = time)
    p = getvar(ncfile, "pressure", timeidx=time)
    geo_500 = smooth2d(interplevel(z, p, 500.), 13, cenweight=1)
    vort_500 = cal_vor(time)[4]
    min_latitude = 145
    max_latitude = 165
    min_longitude = -25
    max_longitude = -12   #era5 is from -5. 
    
    
    
    latitude, longitude = latlon_coords(z_meta)
    
    cart_proj = get_cartopy(z_meta)
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([min_latitude, max_latitude, min_longitude, 
                  max_longitude], crs=crs.PlateCarree())
    max_value = np.nanmax(vort_500)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
  
    levels = np.linspace(5700, 5950, 15)
    vorticity_levels = np.linspace(-2e-3, 2e-3, 21)
    step=20

    plt.contourf(to_np(longitude)[::step, ::step], to_np(latitude)[::step, ::step],
                 to_np(vort_500)[::step, ::step], levels=vorticity_levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label=f'Vorticity ($s^{-1}$)', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    
    # Plot contour lines of vorticity at 200 hPa
    contour =plt.contour(to_np(longitude), to_np(latitude), to_np(geo_500), levels=levels, colors='black',
                transform=crs.PlateCarree())
    plt.clabel(contour, inline=True, fontsize=10, fmt='%.0f m')

    
    # Add other customizations and titles if desired
    plt.title(f'Geopotential Height and Vorticity at 500 hPa {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
   
    save_dir = '/scratch/vj74/ts7941/wrf_geo_500'

    plt.savefig(os.path.join(save_dir,  f"wrf_geo_500 {alphabet[a]}.png"))
    print("plotted geo500",a)
    #plt.show()



def humidity(time,a):
    
    
    
    z = getvar(ncfile, "z", timeidx = time )
    z_meta = getvar(ncfile, "z", timeidx = time)
    p = getvar(ncfile, "pressure", timeidx=time)
    humidity = getvar(ncfile, "rh", timeidx = time)
    humidity_500 = interplevel(humidity, p, 500)
    u = getvar(ncfile, "ua", timeidx = time )
    v = getvar(ncfile, "va", timeidx = time )
    u_500 = interplevel(u, p, 500.)
    v_500 = interplevel(v, p, 500.)
    min_latitude = 145
    max_latitude = 165
    min_longitude = -25
    max_longitude = -12   #era5 is from -5. 
    
    
    
    latitude, longitude = lat, lon
    
    cart_proj = get_cartopy(z_meta)
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([min_latitude, max_latitude, min_longitude, 
                  max_longitude], crs=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    
  
    # Plot shaded contour map of geopotential height
    levels = np.linspace(0, 110, 10)
   

    step = 20
    scale = 15
     
     
    plt.contourf(longitude[::step, ::step], latitude[::step, ::step],
                  humidity_500[::step, ::step], levels=levels, cmap='YlGn',
                  transform=crs.PlateCarree())
    plt.colorbar(label='Relative Humidity (%)', 
             fraction=0.0295, pad=0.04)   
    
    
    
    ax.quiver(to_np(longitude)[::step, ::step], to_np(latitude)[::step, ::step],
              u_500[::step, ::step], v_500[::step, ::step],  transform=crs.PlateCarree(),
              angles='xy', scale_units='xy', scale=scale, width=0.003,
              pivot='middle', edgecolor='k')
    max_value = np.nanmax(humidity_500)
    ax.set_title(f"Max Value: {max_value: .4f}", loc='right')
    
    
    # Add other customizations and titles if desired
    plt.title(f'Relative Humidity at 500 hPa {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
   
    save_dir = '/scratch/vj74/ts7941/wrf_humid'

    plt.savefig(os.path.join(save_dir,  f"wrf_humidity {alphabet[a]}.png"))
    print("plotted humidity", a)
    


def geopotential_200(time,a):
    
   
    
    z = getvar(ncfile, "z", timeidx = time )
    z_meta = getvar(ncfile, "z", timeidx = time)
    p = getvar(ncfile, "pressure", timeidx=time)
    geo_200 = smooth2d(interplevel(z, p, 200.), 13, cenweight=1)
    vort_200 = cal_vor(time)[1]
    min_latitude = 145
    max_latitude = 165
    min_longitude = -25
    max_longitude = -12   #era5 is from -5. 
    
    
    
    latitude, longitude = lat, lon
    
    cart_proj = get_cartopy(z_meta)
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([min_latitude, max_latitude, min_longitude, 
                  max_longitude], crs=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    max_value = np.nanmax(vort_200)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
    vorticity_levels = np.linspace(-2e-3, 2e-3, 21)
    levels = np.linspace(12350, 12550, 15)
    step = 20
    

    plt.contourf(to_np(longitude)[::step, ::step], to_np(latitude)[::step, ::step],
                 to_np(vort_200)[::step, ::step], levels=vorticity_levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label=f'Vorticity ($s^{-1}$)', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    
    # Plot contour lines of vorticity at 200 hPa
    contour = plt.contour(to_np(longitude), to_np(latitude), to_np(geo_200), levels=levels, colors='black',
                transform=crs.PlateCarree())
    
    plt.clabel(contour, inline=True, fontsize=10, fmt='%.0f m')
    # Add other customizations and titles if desired
    plt.title(f'Geopotential Height and Vorticity at 200 hPa {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
   
    save_dir = '/scratch/vj74/ts7941/wrf_geo_200'

    plt.savefig(os.path.join(save_dir,  f"wrf_geo_200 {alphabet[a]}.png"))
    print("plotted geo200",a)
    #plt.show()


def upperwind_and_vort(time,a):
    
    
    u = getvar(ncfile, "ua", timeidx = time )
    v = getvar(ncfile, "va", timeidx = time )
    wspd = smooth2d(getvar(ncfile, "uvmet10", timeidx =time), 3, cenweight=1)
    wspd = np.sqrt(wspd[0]**2 + wspd[1]**2)
    z_meta = getvar(ncfile, "z", timeidx = time)
    
    p = getvar(ncfile, "pressure", timeidx=time)
    
    vort_200 = cal_vor(time)[1]
    u_200 = interplevel(u, p, 200.)
    v_200 = interplevel(v, p, 200.)
        
    
    min_latitude = 145
    max_latitude = 165
    min_longitude = -25
    max_longitude = -12   #era5 is from -5. 
    
    
    
    latitude, longitude = lat, lon
    
    cart_proj = get_cartopy(z_meta)
     
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([min_latitude, max_latitude, min_longitude, 
                  max_longitude], crs=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    max_value = np.nanmax(wspd)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
    # Set a larger scale value to make the arrows longer
    scale = 15
    # Plot shaded contour map of relative vorticity
    levels = np.linspace(-2e-3, 2e-3, 21)

    substitute_levels = np.linspace(220,240,15)
    step = 20
    
    
    plt.contourf(longitude[::step, ::step], latitude[::step, ::step],
                 vort_200[::step, ::step], levels=levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    
    plt.colorbar(label='Relative Vorticity $(s^{-1})$', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    # Plot wind vectors on the map for the first time step at 850 hPa.
    ax.quiver(to_np(longitude)[::step, ::step], to_np(latitude)[::step, ::step],
              u_200[::step, ::step], v_200[::step, ::step],  transform=crs.PlateCarree(),
              angles='xy', scale_units='xy', scale=scale, width=0.003,
              pivot='middle', edgecolor='k')
    
    
    plt.title(f'Wind Vectors and Relative Vorticity at 200 hPa {times[time]}', loc = 'left')
    save_dir = '/scratch/vj74/ts7941/wrf_upperwind'

    plt.savefig(os.path.join(save_dir,  f"wrf_upperwind {alphabet[a]}.png"))
    print("plotted upper",a)
    #plt.show()

def lowerwind_and_vort(time,a):
    
    
    u = getvar(ncfile, "ua", timeidx = time )
    v = getvar(ncfile, "va", timeidx = time )
    z_meta = getvar(ncfile, "z", timeidx = time)
    wspd = smooth2d(getvar(ncfile, "uvmet10", timeidx =time), 3, cenweight=1)
    wspd = np.sqrt(wspd[0]**2 + wspd[1]**2)
    
    p = getvar(ncfile, "pressure", timeidx=time)
    
    u_850 = interplevel(u, p, 850.)
    v_850 = interplevel(v, p, 850.)
    
    
    vort_850 = cal_vor(time)[0]
    
    min_latitude = 145
    max_latitude = 165
    min_longitude = -25
    max_longitude = -12   #era5 is from -5. 
    
    
    
    latitude, longitude = lat, lon
    
    cart_proj = get_cartopy(z_meta)
     
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([min_latitude, max_latitude, min_longitude, 
                  max_longitude], crs=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    max_value = np.nanmax(wspd)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
    # Set a larger scale value to make the arrows longer
    scale = 15
    # Plot shaded contour map of relative vorticity
    levels = np.linspace(-2e-3, 2e-3, 21)
    step = 20
    
    #plt.contourf(longitude[::step, ::step], latitude[::step, ::step], vort_250[::step, ::step], levels=substitute_levels, cmap='coolwarm',
     #            transform=crs.PlateCarree())
    plt.contourf(longitude[::step, ::step], latitude[::step, ::step],
                 vort_850[::step, ::step], levels=levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    
    plt.colorbar(label='Relative Vorticity $(s^{-1})$', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    # Plot wind vectors on the map for the first time step at 850 hPa.
    ax.quiver(to_np(longitude)[::step, ::step], to_np(latitude)[::step, ::step],
              u_850[::step, ::step], v_850[::step, ::step],  transform=crs.PlateCarree(),
              angles='xy', scale_units='xy', scale=scale, width=0.003,
              pivot='middle', edgecolor='k')
    
   
    
    plt.title(f'Wind Vectors and Relative Vorticity at 850 hPa {times[time]}', loc = 'left')

   
    save_dir = '/scratch/vj74/ts7941/wrf_lowerwind'

    plt.savefig(os.path.join(save_dir,  f"wrf_lowerwind {alphabet[a]}.png"))
    print("plotted lower", a)
    
    
    
    
    

def SLP(time, a):
    
    
    
    z = getvar(ncfile, "z", timeidx = time )
    z_meta = getvar(ncfile, "z", timeidx = time)
    
    slp = smooth2d(getvar(ncfile, "slp", timeidx = time), 3, cenweight=1)
    
    
    
    
    p = getvar(ncfile, "pressure", timeidx=time)
    wspd = smooth2d(getvar(ncfile, "uvmet10", timeidx =time), 3, cenweight=1)
    wspd = np.sqrt(wspd[0]**2 + wspd[1]**2)
    
    
    
    min_latitude = 145
    max_latitude = 165
    min_longitude = -25
    max_longitude = -12   #era5 is from -5. 
    
    
    
    latitude, longitude = lat, lon
    
    cart_proj = get_cartopy(z_meta)
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([min_latitude, max_latitude, min_longitude, 
                  max_longitude], crs=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    max_value = np.nanmax(wspd)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
  
    # Plot shaded contour map of geopotential height
    levels = np.linspace(990, 1030, 15)
    wspd_levels = np.linspace(0, 35, 10)

    plt.contourf(to_np(longitude), to_np(latitude), to_np(wspd), levels=wspd_levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label=f'10-m Wind Speed ($ms^{-1}$)', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    
    # Plot contour lines of vorticity at 200 hPa
    contour = plt.contour(to_np(longitude), to_np(latitude), to_np(slp), levels=levels, colors='black',
                transform=crs.PlateCarree())
    plt.clabel(contour, inline=True, fontsize=10, fmt='%.0f hPa')

    # Add other customizations and titles if desired
    plt.title(f'Sea Level Pressure with 10-m Wind speed {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
   
    save_dir = '/scratch/vj74/ts7941/wrf_SLP'

    plt.savefig(os.path.join(save_dir,  f"wrf_slp {alphabet[a]}.png"))
    print("plotted slp", a)
    
    
    
    
    
    
def SST(time, a):
    
    
    
    z = getvar(ncfile, "z", timeidx = time )
    z_meta = getvar(ncfile, "z", timeidx = time)
    
    sst = getvar(ncfile, "T2", timeidx = time)/10
    
    
  
    
    
    min_latitude = 145
    max_latitude = 165
    min_longitude = -25
    max_longitude = -12   #era5 is from -5. 
    
    
    
    latitude, longitude = lat, lon
    
    cart_proj = get_cartopy(z_meta)
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([min_latitude, max_latitude, min_longitude, 
                  max_longitude], crs=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    max_value = np.nanmax(sst)
    ax.set_title(f"Max Value: {max_value:.4f}", loc='right')
  
    # Plot shaded contour map of geopotential height
    levels = np.linspace(20 , 35, 10)

    plt.contourf(to_np(longitude), to_np(latitude), to_np(sst), levels=levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label=f'Sea Surface Temperature (celsius)',
                 fraction=0.0295, pad=0.04)
    
    # Add other customizations and titles if desired
    plt.title(f'Sea Surface Temperature {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
   
    #save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/geo_500'

    save_dir = '/scratch/vj74/ts7941/wrf_SST'

    plt.savefig(os.path.join(save_dir,  f"wrf_sst {alphabet[a]}.png"))
    
    #plt.savefig(os.path.join(save_dir, f"geo500 {alphabet[time]}.png"))
    print("plotted sst",a)
    
    
    
    
    
    
    
   

def wind_shear(time, a):
    
    
    u = getvar(ncfile, "ua", timeidx = time )
    v = getvar(ncfile, "va", timeidx = time )
    z_meta = getvar(ncfile, "z", timeidx = time)
    
    p = getvar(ncfile, "pressure", timeidx=time)
    
    vorticity = getvar(ncfile, "temp", timeidx = time)
    u_850 = interplevel(u, p, 850.)
    v_850 = interplevel(v, p, 850.)
    u_500 = interplevel(u, p, 500.)
    v_500 = interplevel(v, p, 500.)
    
    
    wind_shear = np.sqrt((u_500-u_850)**2 + (v_500-v_850)**2)

    
    
    min_latitude = 145
    max_latitude = 165
    min_longitude = -25
    max_longitude = -12   #era5 is from -5. 
    
    
    
    latitude, longitude = lat, lon
    
    cart_proj = get_cartopy(z_meta)
     
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_extent([min_latitude, max_latitude, min_longitude, 
                  max_longitude], crs=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    max_value = np.nanmax(wind_shear)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
    # Plot shaded contour map of relative vorticity
    shear_levels = np.linspace(0, 40, 10)
    
    plt.contourf(longitude, latitude,
                 wind_shear, levels=shear_levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    
    plt.colorbar(label='Wind Shear $(ms^{-1})$', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    
    
    plt.title(f'Wind Shear {times[time]}', loc = 'left')

    save_dir = '/scratch/vj74/ts7941/windshear'

    plt.savefig(os.path.join(save_dir,  f"wrf_windshear {alphabet[a]}.png"))
    
    print("plotted wind shear", a)
    
    #plt.show() 
    
    
    

