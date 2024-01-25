# created Jan 2024 by G Oldford
# purpose: test for trend at Nanoose, model and observ.
# to-do:
#    - seasonal and non-seasonal test
#    - evaluate the whole Nanoose stn dataset
#    - try averaging to annually and monthly
#    - different depths

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import glob
import cartopy
from cartopy import crs, feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.patches import Rectangle
import cmocean as cm
import os
import time
import datetime as dt

from scipy import interpolate
from scipy.stats import kendalltau
from scipy.stats import linregress
import pymannkendall as mk

import pandas as pd
import mannkendall as mk
import netCDF4 as nc

print(np.sum(1))

var = "temp" # temp or salt
plot_which = 'mk' #mk or lr

proj_p = 'C:/Users/Greig/Documents/github/NEMO-Salish-Sea-1500/'
data_p = ''

modelruns_info = {'SalishSea1500-RUN203': {'path': '../../climatol_intermediate_files/',
                                           'shortcode': 'RUN203',
                                           'label': 'HOTSS-1.1'},
                  'SalishSea1500-RUN216': {'path': '../../climatol_intermediate_files/',
                                           'shortcode': 'RUN216',
                                           'label': 'HOTSS-1.2'}
                  }

depth_min = 4.5
depth_max = 400
mvg_avg = 1
n_iterations = 1000
tyr = np.arange(1980,2019, 1/24)
alpha = 0.05
#
# def reshape_seasonal(data_in):
#
#
#


# get depths, mask, e3t_0 at Nanoose
meshm_p = 'data/mesh mask/'
meshm_f = 'mesh_mask_20210406.nc'
with nc.Dataset(os.path.join(proj_p, meshm_p, meshm_f), 'r') as mesh:
    tmask=mesh.variables['tmask'][:]
    e3t0=mesh.variables['e3t_0'][:]
    gdept_0=mesh.variables['gdept_0'][:]
gdept_0 = gdept_0[0,:,175,75]
e3t0 = e3t0[0,:,175,75]
tmask = tmask[0,:,175,75]

# find closest idx in depth to the min / max depths
diffs_min = np.abs(gdept_0 - depth_min)
min_idx = np.argmin(diffs_min)
diffs_max = np.abs(gdept_0 - depth_max)
max_idx = np.argmin(diffs_max)

e3t0_2 = e3t0[min_idx:max_idx+1]
gdept_0 = gdept_0[min_idx:max_idx+1]

# plot size scalar
fact = 1.4
fig, axs = plt.subplots(3,1, figsize=(5*fact, 5*fact),
                        facecolor='w', edgecolor='k',
                        constrained_layout=True)
axs = axs.ravel()
m = 0



# print(modelruns_info.keys())
# for model in modelruns_info.keys():
#     print('')
#     print('##### MODEL ', model, '#####')
#     mod_run = modelruns_info[model]['shortcode']
#     mod_lbl = modelruns_info[model]['label']
#     path1 = modelruns_info[model]['path']
#
#     obs_anom_f = mod_run + 'obs_anom_' + var + '_1980-2018.nc'
#     mod_anom_f = mod_run + 'mod_anom_' + var + '_1980-2018.nc'
#     obs_anom = xr.open_dataset(os.path.join(path1,obs_anom_f))
#     mod_anom = xr.open_dataset(os.path.join(path1, mod_anom_f))
#
#     obs = obs_anom.isel(deptht = obs_anom.deptht >= depth_min)
#     mod = mod_anom.isel(deptht = mod_anom.deptht >= depth_min)
#     obs = obs.isel(deptht = obs.deptht <= depth_max)
#     mod = mod.isel(deptht = mod.deptht <= depth_max)
#
#     # added this 2024-16
#     obs_dt = obs['time_counter']
#     mod_dt = obs['time_counter']
#
#     # Creating a numpy array with dtype of Python datetime
#     years = obs_dt['time_counter'].dt.year
#     months = obs_dt['time_counter'].dt.month
#     days = obs_dt['time_counter'].dt.day
#     obs_dt_py = np.array([dt.datetime(year, month, day) for year, month, day in zip(years.values, months.values, days.values)])
#
#     if var == 'temp':
#         obs = obs['temperature']
#         mod = mod['temperature']
#     else:
#         obs = obs['salinity']
#         mod = mod['salinity']
#
#     ##################################################
#     # calculate depth integrated mean
#
#     # replace values in e3t where obs are nans with nans
#     nan_indices = np.isnan(obs)
#     # add dimensions to e3t to match shape of data
#     obs_shape = obs.shape  # should be 34, 936
#     e3t0_3 = e3t0_2[:, np.newaxis]
#     e3t0_3 = np.repeat(e3t0_3, obs_shape[1], axis=1)
#     e3t0_3[nan_indices] = np.nan
#
#     # calculate weighted average using e3t_0
#     sum_e3t = np.nansum(e3t0_3, axis=0)
#     e3t0_weights = np.where(sum_e3t != 0, e3t0_3 / sum_e3t[np.newaxis, :], 0)
#     #     crosscheck_weights = np.nansum(e3t0_weights, axis=0)
#     #     print(crosscheck_weights) # should sum to one or be zeros (zeros replaced with nan later)
#     weighted_vals = obs * e3t0_weights
#     obs_avg = np.sum(weighted_vals, axis=0)  # sum instead of nansum keeps it as xarray dataarray
#     obs_avg = obs_avg.where(obs_avg != 0, np.nan)
#     weighted_vals = mod * e3t0_weights
#     mod_avg = np.sum(weighted_vals, axis=0)
#     mod_avg = mod_avg.where(mod_avg != 0, np.nan)
#     ##################################################
#
#     if mvg_avg > 1:
#         modnonan = mod_avg[~np.isnan(mod_avg)]
#         obsnonan = obs_avg[~np.isnan(obs_avg)]
#         mod_avg_mvg = np.convolve(modnonan.data, np.ones(mvg_avg)/mvg_avg, mode='valid')
#         obs_avg_mvg = np.convolve(obsnonan.data, np.ones(mvg_avg)/mvg_avg, mode='valid')
#         mod_avg = mod_avg_mvg
#         obs_avg = obs_avg_mvg
#
#     # run mk_temp_aggr from mannkendall package (not pymannkendall
#     # https://github.com/mannkendall/Python/blob/1dd20639898cfd767c274e8cb75bd6a8b70e3390/src/mannkendall/mk_main.py#L119C5-L119C17
#     # input list of 1-D ndarray of date.datetime with each array defining season
#
#     # seasons - mk requires 1-D list of ndarrays of length = seasons
#     # Grouping by month and creating a list of 1-D ndarrays
#     # seasonal_dates = []
#     # for _, group in obs_dt.groupby('time_counter.month'):
#     #     seasonal_dates.append(np.asarray(group))
#
#     # Define the number of seasons per year
#     seasons_per_year = 24
#     # Grouping by season and creating a list of 1-D ndarrays
#     seasonal_dates = []
#     obs_avg_multi = []
#     mod_avg_multi = []
#     # below does not work b/c of leap years affecting biweekly dates
#     # for seasondt, group in obs_dt.groupby("time_counter.dayofyear"):
#     #     print(group)
#     # # Extracting the year from the datetime array
#     # years = obs_dt['time_counter'].dt.year
#     #
#
#     # pd_timestamp = obs_dt['time_counter']
#     # python_datetime = pd_timestamp.to_pydatetime()
#     # seasons_dt = np.array([python_datetime])
#
#     years = obs_dt['time_counter'].dt.year
#     months = obs_dt['time_counter'].dt.month
#     days = obs_dt['time_counter'].dt.day
#
#
#     # Creating a numpy array with dtype of Python datetime
#     seasons_dt = np.array([dt.datetime(year, month, day) for year, month, day in zip(years.values, months.values, days.values)])
#
#     # this creates biweekly data
#     for season in range(1, seasons_per_year + 1):
#         indices = np.where(
#             (obs_dt['time_counter'].dt.dayofyear >= (season * 15) - 10) &
#             (obs_dt['time_counter'].dt.dayofyear < ((season + 1) * 15) - 10))[0]
#         seasonal_dates.append(np.asarray(seasons_dt[indices]))
#         obs_avg_multi.append(np.asarray(obs_avg[indices]))
#         mod_avg_multi.append(np.asarray(mod_avg[indices]))
#
#     out = mk.mk_temp_aggr(seasonal_dates, obs_avg_multi, 0.00001)
#     # Print the results
#     for n in range(seasons_per_year):
#         print('Observations')
#         print('Season {ind}:'.format(ind=n + 1), out[n])
#
#     out = mk.mk_temp_aggr(seasonal_dates, mod_avg_multi, 0.00001)
#
#     # Print the results
#     for n in range(seasons_per_year):
#         print('Model')
#         print('Season {ind}:'.format(ind=n + 1), out[n])
#
#     out = mk.mk_temp_aggr(obs_dt_py, np.asarray(obs_avg), 0.00001)
#     print(out)
#     out = mk.mk_temp_aggr(obs_dt_py, np.asarray(mod_avg), 0.00001)
#     print(out)
#
#     print("success")

# load data from whole
obs_mc07_f = 'Nanoose_obs_anom_temp_1970-2005.nc'
#obs_mc07_f = 'obs_temperature_1970-2005-bimonthly_timeseries_mc07.nc'
obs_mc07 = xr.open_dataset(os.path.join('../../climatol_intermediate_files/',obs_mc07_f))
obs = obs_mc07.isel(deptht = obs_mc07.deptht >= depth_min)
obs_mc07 = obs.isel(deptht = obs.deptht <= depth_max)

# calculate depth integrated mean
# replace values in e3t where obs are nans with nans
obs_mc07 = obs_mc07['temperature']
nan_indices = np.isnan(obs_mc07)
# add dimensions to e3t to match shape of data

obs_shape = obs_mc07.shape  # should be 34, 936
e3t0_3 = e3t0_2[:, np.newaxis]
e3t0_3 = np.repeat(e3t0_3, obs_shape[1], axis=1)
e3t0_3[nan_indices] = np.nan

# calculate weighted average using e3t_0
sum_e3t = np.nansum(e3t0_3, axis=0)
e3t0_weights = np.where(sum_e3t != 0, e3t0_3 / sum_e3t[np.newaxis, :], 0)
#     crosscheck_weights = np.nansum(e3t0_weights, axis=0)
#     print(crosscheck_weights) # should sum to one or be zeros (zeros replaced with nan later)
weighted_vals = obs_mc07 * e3t0_weights
obs_mc07_avg = np.sum(weighted_vals, axis=0)  # sum instead of nansum keeps it as xarray dataarray
obs_mc07_avg = obs_mc07_avg.where(obs_mc07_avg != 0, np.nan)

obs_dt = obs_mc07['time_counter']
years = obs_dt['time_counter'].dt.year
months = obs_dt['time_counter'].dt.month
days = obs_dt['time_counter'].dt.day
obs_dt_py = np.array([dt.datetime(year, month, day) for year, month, day in zip(years.values, months.values, days.values)])
#

out = mk.mk_temp_aggr(obs_dt_py, np.asarray(obs_mc07_avg), 0.00001)
print('Non-Seasonal 3PW result')
print(out)

# Linear regression
# linear regression needs numbers not dates for time
# convert dates to decimal year
tyr = np.arange(1970,2006,(1/24))
obs_mc07_avg_nonan = obs_mc07_avg[~np.isnan(obs_mc07_avg)]
tyr_nonan = tyr[~np.isnan(obs_mc07_avg)]
A = np.vstack([tyr_nonan, np.ones(len(tyr_nonan))]).T
m, c = np.linalg.lstsq(A, np.asarray(obs_mc07_avg_nonan), rcond=None)[0]
print("Linear Regression result:")
print('Slope: ' + str(m))

# seasonal
# this creates 1-D list of ndarray dates as required by mannkendall
seasons_per_year = 24
seasonal_dates = []
obs_mc07_avg_multi = []
for season in range(1, seasons_per_year + 1):
    indices = np.where(
        (obs_dt.dt.dayofyear >= (season * 15) - 10) &
        (obs_dt.dt.dayofyear < ((season + 1) * 15) - 10))[0]
    seasonal_dates.append(np.asarray(obs_dt_py[indices]))
    obs_mc07_avg_multi.append(np.asarray(obs_mc07_avg[indices]))
print(len(seasonal_dates))
print(len(obs_mc07_avg_multi))
out = mk.mk_temp_aggr(seasonal_dates, obs_mc07_avg_multi, 0.00001)
print('Seasonal 3PW Result: ')
for n in range(seasons_per_year):
    print('Season {ind}:'.format(ind=n + 1), out[n])



# aggregate to annual
obs_annual_avg_mc07 = obs_mc07_avg.groupby("time_counter.year").mean(dim="time_counter")
obs_dt_annual_mc07 = obs_annual_avg_mc07['year']
years = obs_dt_annual_mc07['year']
obs_dt_annual_mc07_py = np.array([dt.datetime(year,1,1) for year in years.values])
out = mk.mk_temp_aggr(obs_dt_annual_mc07_py, np.asarray(obs_annual_avg_mc07), 0.00001)

print('1970 - 2005, annually averaged, Nanoose Obs, Depth 4 - 400 m')
print('Non-Seasonal 3PW result')
print(out)



# obs_mc07_avg
# for yr in range(1970,2006)
#     filt =
#

    # tperiod = ds_mc07.timeperiod.values
    # timestamp = ds_mc07.time_counter.values
    # dates_obs = pd.DatetimeIndex(timestamp)
    # obs_yrs = np.array(dates_obs.year)
    #
    # salt = ds_mc07.salinity.values
    # temp = ds_mc07.temperature.values
    #
    # salt_bimonth = np.zeros([864,40]); salt_bimonth[:] = np.nan
    # temp_bimonth = np.zeros([864,40]); temp_bimonth[:] = np.nan
    # this_timeperiod = np.zeros([864])

    # ind = 0
    # temp_bimonth = np.zeros([864, 40]);
    # temp_bimonth[:] = np.nan
    # for yr in range(1970,2006):
    #     filt = ((obs_yrs == yr) & (tperiod == tp))
    #     temp_bimonth[ind,:] = np.nanmean((temp[filt,:]), axis = 0)
    #     this_timeperiod[ind] = tp
    #     ind = ind+1