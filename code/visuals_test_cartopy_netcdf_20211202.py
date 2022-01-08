# created Dec 1 2021 by GO
# - issues with applying projections 
# - following http://tech.weatherforce.org/blog/ecmwf-data-animation/index.html
# - you can't add cartopy layers for background etc while not declaring projection, but declaring projection doesn't work with my lat lon data


import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.animation as animation
import netCDF4 as nc
import numpy as np
from netCDF4 import Dataset


### ANIMATE CODE FOR NEMO (using xarray)
def make_figure_NEMO(xlim, ylim):
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    return fig, ax
    
def draw_NEMO(frame, add_colorbar):
    field_at_time = field[frame]
    contour = field_at_time.plot(ax=ax, add_colorbar=add_colorbar, vmin=min_value, vmax=max_value)
    #title = u"%s — %s" % ("field name here", str(ds.time_counter[frame].values)[:19])
    #ax.set_title(title)
    return contour
    
def animate_NEMO(frame):
    return draw_NEMO(frame, add_colorbar=False)
    
    
def init_NEMO():
    return draw_NEMO(0, add_colorbar=True)



### ANIMATE CODE FOR ERA (using netcdf4)
def make_figure_ERA():
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(1, 1, 1)
    
    return fig, ax

def init_ERA():
    return draw_ERA(0, add_colorbar=True)
    
def draw_ERA(frame, add_colorbar):
    field_at_time = field[frame]
    ax.pcolormesh(field_at_time)

    #title = u"%s — %s" % ("field name here", str(ds.time_counter[frame].values)[:19])
    #ax.set_title(title)
    return ax

def animate_ERA(frame):
    return draw_ERA(frame, add_colorbar=False)


# original code
def draw_orig(frame, add_colorbar):
    field_at_time = ds.mldkz5[frame]
    contour = field_at_time.plot(ax=ax, transform=ccrs.PlateCarree(),
                        add_colorbar=add_colorbar, vmin=min_value, vmax=max_value)
    title = u"%s — %s" % ("field name here", str(ds.time_counter[frame].values)[:19])
    ax.set_title(title)
    return contour

def make_figure():
    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # generate a basemap with country borders, oceans and coastlines
    ax.add_feature(cfeat.LAND)
    ax.add_feature(cfeat.OCEAN)
    ax.add_feature(cfeat.COASTLINE)
    ax.add_feature(cfeat.BORDERS, linestyle='dotted')
    return fig, ax

make_figure();


#############################
######## NEMO Outputs #######
#############################
nc_file_path = "..//data//temp//"
nc_file = "SalishSea1500-RUN201_MonthlyMeans_grid_T_2D_1981.nc"

ds = xr.open_mfdataset(nc_file_path + nc_file) #open_mfdataset and open_dataset are different

field = ds.mldkz5 # turbocline

max_lat = field.nav_lat.values.max()
min_lat = field.nav_lat.values.min()
max_lon = field.nav_lon.values.max()
min_lon = field.nav_lon.values.min()
print(max_lon)

print("turbocline data retrieved")
print(field.shape)

xlim = (min_lon, max_lon)
ylim = (min_lat, max_lat)
_, ax = make_figure_NEMO(xlim, ylim)

#field_at_time.plot(ax=ax, transform=ccrs.PlateCarree(), vmin=min_value, vmax=max_value)
field_at_time = field[3]
print(field_at_time.shape)
field_at_time.plot(ax = ax)

plt.savefig('nemo_visual_test_mldkz5.png')
plt.close()

ax.cla()
print("printed NEMO turbocline (mldkz5) image")

############ animation ###############
#frames = ds.time.size 
##frames = ds.time_counter.size        # Number of frames
#print("# of frames" + str(frames))

#ax.cla()
#fig, ax = make_figure_GO(xlim, ylim)

#print("generating animation")
#ani = animation.FuncAnimation(fig, animate, frames, interval=0.01, blit=False,
#                              init_func=init, repeat=False)
#ani._init_draw()

#print("saving animation")
#fps_go=2
#ani.save('test.mp4', writer=animation.FFMpegWriter(fps=fps_go))

#print("success")
#plt.close(fig)


##################################
######## REANALYSIS DATA #########
##################################

# annual file hourly data trimmed to Salish Sea 
#nc_file_path = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//29. Oceanographic Atmospheric//ECMWF ERA5//adjusted//"
#nc_file="/ERA5_SalishSea_fixed_1981.nc"

# monthly means  
nc_file_path = r"C:/Users/Greig/Sync/For Joe/"
nc_file="ERA5_NEMOgrid_light_monthly_1981.nc"
#nc_file="ERA5_NEMOgrid_light_daily_2019.nc"

with nc.Dataset(nc_file_path + nc_file) as ncf:
    nlat = ncf.variables['latitude'][0,...]
    nlon = ncf.variables['longitude'][0,...]
    msdwswrf = ncf.variables['msdwswrf'][:,...] # incoming light
    time1 = ncf.variables['time'][:,...]

print("Light data retrieved from reanalysis product")
print(msdwswrf.shape)

field = msdwswrf

max_lat = nlat.max()
min_lat = nlat.min()
max_lon = nlon.max()
min_lon = nlon.min()
print(min_lon)

xlim = (min_lon, max_lon)
ylim = (min_lat, max_lat)
_, ax = make_figure_ERA()

#field_at_time.plot(ax=ax, transform=ccrs.PlateCarree(), vmin=min_value, vmax=max_value)
#print(np.squeeze(field[0]).shape)
field_at_time = field[0]
ax.pcolormesh(field_at_time)
plt.savefig('nemo_visual_test_msdwswrf.png')
plt.close()
ax.cla()
print("printed ERA shortwave (msdwswrf) image")

#======================================
# ============= animation =============
frames = len(time1)

ax.cla()
fig, ax = make_figure_ERA()

print("generating animation")
print("length: " + str(frames))
ani = animation.FuncAnimation(fig, animate_ERA, frames, interval=0.01, blit=False,
                              init_func=init_ERA, repeat=False)
#ani._init_draw()

print("saving animation")
fps_go=1
ani.save('era_animation_test.mp4', writer=animation.FFMpegWriter(fps=fps_go))

print("success")
plt.close(fig)
#======================================

##########################################
########### ECOSPACE ASC DATA ############
##########################################

