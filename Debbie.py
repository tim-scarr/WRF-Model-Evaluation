# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 15:35:14 2023

@author: timsc
"""

from netCDF4 import Dataset
from netCDF4 import num2date
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.pyplot import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import cartopy.feature as cfeature
import datetime
from mpl_toolkits import basemap
import pandas
from moviepy.editor import ImageSequenceClip
import os
import numpy as np
import astropy.units as u

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, CoordPair)
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import string
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
import requests
from bs4 import BeautifulSoup
import pandas as pd


ncfile = Dataset('era5_2.nc')




ncfile2 = Dataset('era5_2.nc')

timeun = ncfile.variables['time'].units
time_values = ncfile.variables['time'][:]
# Convert time_values to datetime objects using correct unit string format
date_times = num2date(time_values, units='hours since 1900-01-01 00:00:00.0', calendar='gregorian')
times = [dt.strftime('%Y-%m-%d %H:%M') for dt in date_times]

era_times = [i for i in range(24)]


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)



ncfile3 = Dataset('era5_single_levels_2.nc')

min_latitude = -25
max_latitude = -12
min_longitude = 145
max_longitude = 165.2

def upper_windmap():
    # Assuming you have already loaded the data, you can access U and V components like this:
    u_wind = ncfile.variables['u'][0, 3, :]  # U component at the first time step and 850 hPa (index 3).
    v_wind = ncfile.variables['v'][0, 3, :]  # V component at the first time step and 850 hPa (index 3).
    
    # Assuming you have also loaded the latitude and longitude data.
    latitude = ncfile.variables['latitude'][:]
    longitude = ncfile.variables['longitude'][:]
    
    # Calculate wind speed from U and V components
    wind_speed = np.sqrt(u_wind**2 + v_wind**2)
    
    # Define latitude and longitude boundaries for the plot (same as the boundaries used for downloading ERA5 data)
   
    
    # Subsample the data by selecting every nth data point (increase n for fewer arrows)
    subsample_n = 4
    u_subsampled = u_wind[::subsample_n, ::subsample_n]
    v_subsampled = v_wind[::subsample_n, ::subsample_n]
    wind_speed_subsampled = wind_speed[::subsample_n, ::subsample_n]
    plot_latitude = latitude[::subsample_n]
    plot_longitude = longitude[::subsample_n]
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    # Convert latitude and longitude to map coordinates.
    x, y = np.meshgrid(plot_longitude, plot_latitude)
    
    # Set a larger scale value to make the arrows longer
    scale = 10
    
    # Plot wind vectors on the map for the first time step at 850 hPa.
    ax.quiver(x, y, u_subsampled, v_subsampled,  transform=crs.PlateCarree(),
              angles='xy', scale_units='xy', scale=scale, width=0.003,
              pivot='middle', edgecolor='k')
    
    plt.title('ERA5 Wind Vectors at 850 hPa (00:00 21/03/2017)')










def create_movie_from_images(image_folder, output_file, fps=10):
    # Get the list of PNG files in the specified folder
    image_files = [f for f in os.listdir(image_folder) if f.endswith('.png')]
    image_files.sort()  # Sort the files to maintain order
    
    # Create a list of image file paths
    image_paths = [os.path.join(image_folder, image_file) for image_file in image_files]
    
    # Create the movie clip using ImageSequenceClip
    clip = ImageSequenceClip(image_paths, fps=fps)
    
    # Save the movie as an MP4 file
    clip.write_videofile(output_file)








def lowerwind_and_vort(time):
    u_wind = ncfile.variables['u'][time, 2, :]  # U component at the first time step and 850 hPa (index 3).
    v_wind = ncfile.variables['v'][time, 2, :]  # V component at the first time step and 850 hPa (index 3).
    
    latitude = ncfile.variables['latitude'][:]
    longitude = ncfile.variables['longitude'][:]
    
    # Calculate wind speed from U and V components
    wind_speed = np.sqrt(u_wind**2 + v_wind**2)
    
    
    # Subsample the data by selecting every nth data point (increase n for fewer arrows)
    subsample_n = 5
    u_subsampled = u_wind[::subsample_n, ::subsample_n]
    v_subsampled = v_wind[::subsample_n, ::subsample_n]
    wind_speed_subsampled = wind_speed[::subsample_n, ::subsample_n]
    plot_latitude = latitude[::subsample_n]
    plot_longitude = longitude[::subsample_n]
    
    # Access relative vorticity data from the netCDF file
    relative_vorticity = ncfile.variables['vo'][time, 2, ::subsample_n, ::subsample_n]
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    # Convert latitude and longitude to map coordinates.
    x, y = np.meshgrid(plot_longitude, plot_latitude)
    
    # Set a larger scale value to make the arrows longer
    scale = 15
    # Plot shaded contour map of relative vorticity
    levels = np.linspace(-2e-3, 2e-3, 21)
    plt.contourf(longitude[::subsample_n], latitude[::subsample_n], relative_vorticity, levels=levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    
    plt.colorbar(label='Relative Vorticity $(s^{-1})$', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    # Plot wind vectors on the map for the first time step at 850 hPa.
    ax.quiver(x, y, u_subsampled, v_subsampled,  transform=crs.PlateCarree(),
              angles='xy', scale_units='xy', scale=scale, width=0.003,
              pivot='middle', edgecolor='k')
    
    max_value = np.nanmax(ncfile.variables['vo'][time, 2, ::subsample_n, ::subsample_n])
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
    
    plt.title(f'ERA5 Wind Vectors and Relative Vorticity at 850 hPa {times[time]}', loc = 'left')
    save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/lower_wind_and_vort'

   
    plt.savefig(os.path.join(save_dir, f"lower_wind {alphabet[time]}.png"))
    #plt.show()
    
def upperwind_and_vort(time):
    u_wind = ncfile.variables['u'][time, 0, :]  # U component at the first time step and 850 hPa (index 3).
    v_wind = ncfile.variables['v'][time, 0, :]  # V component at the first time step and 850 hPa (index 3).
    
    latitude = ncfile.variables['latitude'][:]
    longitude = ncfile.variables['longitude'][:]
    
    # Subsample the data by selecting every nth data point (increase n for fewer arrows)
    subsample_n = 5
    u_subsampled = u_wind[::subsample_n, ::subsample_n]
    v_subsampled = v_wind[::subsample_n, ::subsample_n]
    plot_latitude = latitude[::subsample_n]
    plot_longitude = longitude[::subsample_n]
    
    # Access relative vorticity data from the netCDF file
    relative_vorticity = ncfile.variables['vo'][time, 0, ::subsample_n, ::subsample_n]
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    
    # Convert latitude and longitude to map coordinates.
    x, y = np.meshgrid(plot_longitude, plot_latitude)
    
    # Set a larger scale value to make the arrows longer
    scale = 15
    # Plot shaded contour map of relative vorticity
    levels = np.linspace(-2e-3, 2e-3, 21)
    plt.contourf(longitude[::subsample_n], latitude[::subsample_n], relative_vorticity, levels=levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label='Relative Vorticity $(s^{-1})$', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    # Plot wind vectors on the map for the first time step at 850 hPa.
    ax.quiver(x, y, u_subsampled, v_subsampled,  transform=crs.PlateCarree(),
              angles='xy', scale_units='xy', scale=scale, width=0.003,
              pivot='middle', edgecolor='k')
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    max_value = np.nanmax(relative_vorticity)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
    
    plt.title(f'ERA5 Wind Vectors and Relative Vorticity at 200 hPa {times[time]}', loc = 'left')
    save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/upperwind_and_vort'

   
    plt.savefig(os.path.join(save_dir, f"upper_wind {alphabet[time]}.png"))
    #plt.show()    
    
    

def midwind_and_vort():
    u_wind = ncfile.variables['u'][22, 1, :]  # U component at the first time step and 850 hPa (index 3).
    v_wind = ncfile.variables['v'][22, 1, :]  # V component at the first time step and 850 hPa (index 3).
    
    latitude = ncfile.variables['latitude'][:]
    longitude = ncfile.variables['longitude'][:]
    
    # Calculate wind speed from U and V components
    wind_speed = np.sqrt(u_wind**2 + v_wind**2)
    
    # Define latitude and longitude boundaries for the plot (same as the boundaries used for downloading ERA5 data)
    
    
    # Subsample the data by selecting every nth data point (increase n for fewer arrows)
    subsample_n = 5
    u_subsampled = u_wind[::subsample_n, ::subsample_n]
    v_subsampled = v_wind[::subsample_n, ::subsample_n]
    wind_speed_subsampled = wind_speed[::subsample_n, ::subsample_n]
    plot_latitude = latitude[::subsample_n]
    plot_longitude = longitude[::subsample_n]
    
    # Access relative vorticity data from the netCDF file
    relative_vorticity = ncfile.variables['vo'][22, 1, ::subsample_n, ::subsample_n]
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    # Convert latitude and longitude to map coordinates.
    x, y = np.meshgrid(plot_longitude, plot_latitude)
    
    # Set a larger scale value to make the arrows longer
    scale = 15
    # Plot shaded contour map of relative vorticity
    levels = np.linspace(-4e-4, 1.2e-4, 15)
    plt.contourf(longitude[::subsample_n], latitude[::subsample_n], relative_vorticity, levels=levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label='Relative Vorticity $(s^{-1})$')
    
    # Plot wind vectors on the map for the first time step at 850 hPa.
    ax.quiver(x, y, u_subsampled, v_subsampled,  transform=crs.PlateCarree(),
              angles='xy', scale_units='xy', scale=scale, width=0.003,
              pivot='middle', edgecolor='k')
    
   
    
    plt.title('ERA5 Wind Vectors and Relative Vorticity at 200 hPa (12:00 26/03/2017)')
    #plt.savefig("upperwind_and_vort 12 .png")
    plt.show()





def humidity(time):
    u_wind = ncfile.variables['u'][time, 1, :]  # U component at the first time step and 850 hPa (index 3).
    v_wind = ncfile.variables['v'][time, 1, :]  # V component at the first time step and 850 hPa (index 3).
    
    # Assuming you have also loaded the latitude and longitude data.
    latitude = ncfile.variables['latitude'][:]
    longitude = ncfile.variables['longitude'][:]
    
    # Calculate wind speed from U and V components
    wind_speed = np.sqrt(u_wind**2 + v_wind**2)
    
   
    
    # Subsample the data by selecting every nth data point (increase n for fewer arrows)
    subsample_n = 5
    u_subsampled = u_wind[::subsample_n, ::subsample_n]
    v_subsampled = v_wind[::subsample_n, ::subsample_n]
    plot_latitude = latitude[::subsample_n]
    plot_longitude = longitude[::subsample_n]
    
    # Access relative vorticity data from the netCDF file
    relative_humidity = ncfile.variables['r'][time, 1, ::subsample_n, ::subsample_n]
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    
    # Convert latitude and longitude to map coordinates.
    x, y = np.meshgrid(plot_longitude, plot_latitude)
    
    # Set a larger scale value to make the arrows longer
    scale = 15
    # Plot shaded contour map of relative vorticity
    levels = np.linspace(0, 110, 11)
    plt.contourf(longitude[::subsample_n], latitude[::subsample_n], relative_humidity, levels=levels, cmap='YlGn',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label='Relative Humidity (%)',
                 fraction=0.0295, pad=0.04)
    
    # Plot wind vectors on the map for the first time step at 850 hPa.
    ax.quiver(x, y, u_subsampled, v_subsampled,  transform=crs.PlateCarree(),
              angles='xy', scale_units='xy', scale=scale, width=0.003,
              pivot='middle', edgecolor='k')
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    max_value = np.nanmax(relative_humidity)
    ax.set_title(f"Max Value: {max_value: .4f}", loc='right')
    
    plt.title(f'ERA5 Wind Vectors and Relative Humidity at 500 hPa {times[time]}', loc = 'left')
    save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/mid_rel_humid'

   
    plt.savefig(os.path.join(save_dir, f"midrel {alphabet[time]}.png"))

   # plt.show()



def geopotential_500(time):
    geo = ncfile2.variables['z'][time, 1, :]/(9.80665*u.m/u.s**2)  # geopotential height
    vorticity_500hPa = ncfile2.variables['vo'][time, 1, :, :]
    
    
    latitude = ncfile2.variables['latitude'][:]
    longitude = ncfile2.variables['longitude'][:]
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    max_value = np.nanmax(vorticity_500hPa)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
  
    # Plot shaded contour map of geopotential height
    levels = np.linspace(5700, 5950, 15)
    vorticity_levels = np.linspace(-2e-3, 2e-3, 21)

    plt.contourf(longitude, latitude, vorticity_500hPa, levels=vorticity_levels, cmap='coolwarm',
                 transform=crs.PlateCarree())

    
    plt.colorbar(label=f'Vorticity ($s^{-1}$)', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    
    # Plot contour lines of vorticity at 200 hPa
    contour = plt.contour(longitude, latitude, geo, levels=levels, colors='black',
                transform=crs.PlateCarree())
    plt.clabel(contour, inline=True, fontsize=10, fmt='%.0f m')

    
    # Add other customizations and titles if desired
    plt.title(f'Geopotential Height and Vorticity at 500 hPa {times[time]}', loc = 'left')
    
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
   
    save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/geo_500'

   
    plt.savefig(os.path.join(save_dir, f"geo500 {alphabet[time]}.png"))
    #plt.show()
    
    
def geopotential_200(time):
    geo = ncfile2.variables['z'][time, 0, :]/(9.80665*u.m/u.s**2)  # geopotential height
    vorticity_500hPa = ncfile2.variables['vo'][time, 0, :, :]
    
    
    latitude = ncfile2.variables['latitude'][:]
    longitude = ncfile2.variables['longitude'][:]
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    max_value = np.nanmax(vorticity_500hPa)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
  
    # Plot shaded contour map of geopotential height
    vorticity_levels = np.linspace(-2e-3, 2e-3, 21)
    levels = np.linspace(12350, 12550, 15)
    plt.contourf(longitude, latitude, vorticity_500hPa, levels=vorticity_levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label=f'Vorticity ($s^{-1}$)', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    
    # Plot contour lines of vorticity at 200 hPa
    contour = plt.contour(longitude, latitude, geo, levels=levels, colors='black',
                transform=crs.PlateCarree())
    
    plt.clabel(contour, inline=True, fontsize=10, fmt='%.0f m')

    
    # Add other customizations and titles if desired
    plt.title(f'Geopotential Height and Vorticity at 200 hPa {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/geo_200'

   
    plt.savefig(os.path.join(save_dir, f"geo200 {alphabet[time]}.png"))
    
    
    
    
    
    
    
def SLP(time):
    
    
    slp = ncfile3.variables['msl'][time, :]/100
    u_wind = ncfile3.variables['u10'][time] 
    v_wind = ncfile3.variables['v10'][time] 
    
    
    wind_speed = np.sqrt(u_wind**2 + v_wind**2)
    
   
    latitude = ncfile3.variables['latitude'][:]
    longitude = ncfile3.variables['longitude'][:]
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    max_value = np.nanmax(wind_speed)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
  
    # Plot shaded contour map of geopotential height
    levels = np.linspace(990, 1030, 15)
    wspd = np.linspace(0, 35, 10)
    plt.contourf(longitude, latitude, wind_speed, levels=wspd, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label='Wind Speed (m/s)', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    
    # Plot contour lines of vorticity at 200 hPa
    contour = plt.contour(longitude, latitude, slp, levels=levels, colors='black',
                transform=crs.PlateCarree())
    
    plt.clabel(contour, inline=True, fontsize=10, fmt='%.0f hPa')
    
    
    
    
    
    
   
    # Add other customizations and titles if desired
    plt.title(f'Sea Level Pressure with 10-m Wind Speed {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/slp'
   
    plt.savefig(os.path.join(save_dir, f"slp {alphabet[time]}.png"))
    
    
    
    
    




alphabet = list(string.ascii_lowercase)






#input_folder = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/wrf_upperwind'
#output_movie_file = 'wrf_upperwind.mp4'
#create_movie_from_images(input_folder, output_movie_file, fps=2)



    
    
    

    
def wind_shear(time):
    u_850 = ncfile2.variables['u'][time, 2, :]
    v_850 = ncfile2.variables['v'][time, 2, :]
    u_500 = ncfile2.variables['u'][time, 1, :]
    v_500 = ncfile2.variables['v'][time, 1, :]
    
    wind_shear = np.sqrt((u_500-u_850)**2 + (v_500-v_850)**2)
    
    
    
    updrafts = ncfile2.variables['vo'][time,2,:]
    
    
    
    
    latitude = ncfile2.variables['latitude'][:]
    longitude = ncfile2.variables['longitude'][:]
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    max_value = np.nanmax(wind_shear)
    ax.set_title(f"Max Value: {max_value:.4e}", loc='right')
    shear_levels = np.linspace(0, 40, 10)
    
  
    # Plot shaded contour map of geopotential height
    levels = np.linspace(-4e-4, 4e-4 , 15)
    plt.contourf(longitude, latitude, wind_shear, levels=shear_levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label='Wind Shear (m/s)', format=ticker.FuncFormatter(fmt),
                 fraction=0.0295, pad=0.04)
    
    
    
    # Plot contour lines of vorticity at 200 hPa
 
    
    # Add other customizations and titles if desired
    plt.title(f'Wind Shear {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    #plt.plot()
   
    save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/Wind Shear'

   
    plt.savefig(os.path.join(save_dir, f"shear {alphabet[time]}.png"))
    
    
   
    
def SST(time):
    
    
    sst = ncfile3.variables['sst'][time, :]-273.15
    
   
    
   
    latitude = ncfile3.variables['latitude'][:]
    longitude = ncfile3.variables['longitude'][:]
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    max_value = np.nanmax(sst)
    ax.set_title(f"Max Value: {max_value:.4f}", loc='right')
  
    # Plot shaded contour map of geopotential height
    levels = np.linspace(20, 35, 10)
    plt.contourf(longitude, latitude, sst, levels=levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label='Sea Surface Temperature (Celsius)',fraction=0.0295, pad=0.04)
    
    
    
   
    # Add other customizations and titles if desired
    plt.title(f'Sea Surface Temperature  at {times[time]}', loc = 'left')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data'
   
    plt.savefig(os.path.join(save_dir, f"SST.png"))   
    
def updrafts(time):
    updrafts = ncfile2.variables['w'][time,2,:]
    
    
    
    
    latitude = ncfile2.variables['latitude'][:]
    longitude = ncfile2.variables['longitude'][:]
    
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    shear_levels = np.linspace(-4, 4, 15)
    
  
    # Plot shaded contour map of geopotential height
    levels = np.linspace(-4e-4, 1.2e-4 , 10)
    plt.contourf(longitude, latitude, updrafts, levels=shear_levels, cmap='coolwarm',
                 transform=crs.PlateCarree())
    
    plt.colorbar(label='Convective Updrafts Pa / s')
    
    
    
    # Plot contour lines of vorticity at 200 hPa
 
    
    # Add other customizations and titles if desired
    plt.title(f'Vorticity')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    plt.plot()
   
    #save_dir = 'C:/Users/timsc/OneDrive - Australian National University/Internship Data/geo_500'

   
    #plt.savefig(os.path.join(save_dir, f"geo500 {alphabet[time]}.png"))
    
    
import html5lib
    
    

url =    "https://ncics.org/ibtracs/index.php?name=v04r00-2017082S14152"
response = requests.get(url)

# Check if the request was successful (status code 200).
if response.status_code == 200:
    # Parse the HTML content of the page using BeautifulSoup.
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find the table by its name attribute.
    #table_name = 'idata'  # Replace with the actual name attribute.
    table = soup.find('table')

    # Use Pandas to read the table into a DataFrame.
    tables = pd.read_html(url)
    desired_table = tables[3]
    # Now, df contains the table data as a DataFrame.
    #print(desired_table)
#else:
    #print(f"Failed to retrieve the webpage. Status code: {response.status_code}")

    
min_latitude = -30
max_latitude = -12
min_longitude = 140
max_longitude = 160.2
    
ibtracs_data = desired_table.drop(0).apply(pd.to_numeric, errors = 'coerce')

lat_deb = ibtracs_data['LAT']
lon_deb = ibtracs_data['LON']
wspd_deb = ibtracs_data["USA WIND"][:30]
time_deb = ibtracs_data["ISO_TIME_________"][:30]


#wrf_track_headers = ["TIME", "LAT", "LON", "PRESSURE", "WSPD" ]


track_wrf = pd.read_csv("TCcenter.csv")

track_wrf  = track_wrf.iloc[24::12]
#track_wrf = track_wrf.iloc[::2


time_wrf = track_wrf["TIME"]
time_wrf_short = []

for original_string in time_wrf:
    # Perform slicing on each string as needed
    sliced_string = original_string[8:16]  # Example: slice from index 7 to 11
    time_wrf_short.append(sliced_string)

wspd_wrf = track_wrf["WSPD"]

pressure_deb = ibtracs_data["USA PRES"][:30]
pressure_wrf = track_wrf['PRESSURE']


track_wrf2 = pd.read_csv("TCcenter_off.csv")

track_wrf2  = track_wrf2.iloc[24::12]
#track_wrf = track_wrf.iloc[::2


wspd_wrf2 = track_wrf2["WSPD"]

pressure_wrf2 = track_wrf2['PRESSURE']


    
def debbie_track(): 
    
    # Create a map using Cartopy.
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.set_xlim([min_longitude, max_longitude])
    ax.set_ylim([min_latitude, max_latitude])
    plt.plot(track_wrf['LON'][3:], track_wrf['LAT'][3:], marker='o', linestyle='-', color='blue', label = "WRF Ocean Off")
    plt.plot(lon_deb, lat_deb, marker='o', linestyle='-', color='green', label = "IBTRACS")
    plt.plot(track_wrf2['LON'][3:], track_wrf2['LAT'][3:], marker='o', linestyle='-', color='red', label = "WRF Ocean On")
    plt.plot(lon_deb.iloc[0], lat_deb.iloc[0], marker='o', linestyle='-', color='pink')
    plt.plot(track_wrf['LON'].iloc[3], track_wrf['LAT'].iloc[3], marker='o', linestyle='-', color='pink')
    plt.plot(track_wrf2['LON'].iloc[3], track_wrf2['LAT'].iloc[3], marker='o', linestyle='-', color='pink')
    ax.grid(True)
    plt.title("Track Comparison", fontsize = 24)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.legend(fontsize = 'x-large')

    plt.savefig("Debbie Track.png")
    plt.show()
    
    
    
def wspd_intensity():
    plt.figure(figsize=(14, 8))

    plt.plot(time_wrf_short, wspd_wrf, label = "WRF Ocean Off" , color = "blue", linewidth =5)
    plt.plot(time_wrf_short, wspd_deb, label = "IBTRACS" , color = "Green", linewidth = 5)
    plt.plot(time_wrf_short, wspd_wrf2, label = "WRF Ocean On" , color = "red", linewidth =5)
    plt.xticks(rotation=75, fontsize = 16)
    plt.yticks(fontsize = 24)
    plt.xlabel("Date", fontsize = 24)
    plt.ylabel(r"Intensity $(kt)$ ", fontsize = 24)
    plt.title("Intensity Comparison of Wind Speed", fontsize = 24)
    plt.legend(fontsize = "xx-large")
    plt.savefig("Wind Speed Analysis.png")
    plt.show()
    
    
    
    
def pres_intensity():
    plt.figure(figsize=(14, 8))

    plt.plot(time_wrf_short, pressure_wrf, label = "WRF Ocean Off" , color = "blue", linewidth =5)
    plt.plot(time_wrf_short, pressure_deb, label = "IBTRACS" , color = "Green", linewidth = 5)
    plt.plot(time_wrf_short, pressure_wrf2, label = "WRF Ocean On" , color = "red", linewidth = 5)

    plt.xticks(rotation=75, fontsize = 16)
    plt.yticks(fontsize=24)
    plt.xlabel("Date", fontsize = 24)
    plt.ylabel(r"Intensity (mb) ", fontsize = 24)
    plt.title("Intensity Comparison of Pressure", fontsize = 24)
    plt.legend(fontsize= "xx-large")
    plt.savefig("Pressure Analysis.png")
    plt.show()    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    