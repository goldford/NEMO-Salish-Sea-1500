# misc tools by GO for the trend analysis
# by GO Jan 2024

import numpy as np
import netCDF4 as nc
import os
from scipy.fft import fft, ifft
import xarray as xr
import pandas as pd

def do_seasonal_avg(dat_davg, time_inc):
    if time_inc == "annual":
        dat_davg = dat_davg.groupby("time_counter.year").mean(dim="time_counter")
        date_range = [np.datetime64(f'{year}-01-01') for year in dat_davg['year'].values]
        dat_davg['year'] = date_range
        dat_davg = dat_davg.rename({'year': 'time_counter'})
        time_nm = 'time_counter'
        dat_ts = 1
    elif time_inc == "monthly":
        dat_davg['time_counter'] = dat_davg['time_counter'].dt.strftime('%Y-%m')
        dat_davg = dat_davg.groupby('time_counter').mean(dim='time_counter')
        date_range = pd.to_datetime(dat_davg['time_counter'], format='%Y-%m')
        dat_davg['time_counter'] = date_range
        dat_ts = 12
        time_nm = 'time_counter'
    elif time_inc == "seasonal":
        # awkward but no better way found
        seasons = {'winter': [12, 1, 2], 'spring': [3, 4, 5], 'summer': [6, 7, 8], 'fall': [9, 10, 11]}
        dat_davg['time_counter_seas'] = dat_davg['time_counter'].values[0]
        start_yr = dat_davg['time_counter'].dt.year.values[0]
        end_yr = dat_davg['time_counter'].dt.year.values[-1]
        for yr in range(start_yr, end_yr + 1, 1):
            for season_name, season_months in seasons.items():
                date_value = np.datetime64(f'{yr}-{season_months[1]:02}-01')

                dat_davg['time_counter_seas'] = xr.where(
                    (dat_davg['time_counter.month'].isin(season_months) & (dat_davg['time_counter.year'] == yr)),
                    date_value,
                    dat_davg['time_counter_seas'])

        dat_davg = dat_davg.groupby('time_counter_seas').mean(dim='time_counter')

        # second loop to add a 'season' label for later
        dat_davg['season'] = 'other'  # initialise
        for yr in range(start_yr, end_yr, 1):
            for season_name, season_months in seasons.items():
                dat_davg['season'] = xr.where(
                    dat_davg['time_counter_seas.month'].isin(season_months),
                    season_name,
                    dat_davg['season'])

        dat_ts = 4
        time_nm = 'time_counter_seas'
    else:
        print("data are assumed biweekly")
        dat_ts = 24
        time_nm = 'time_counter'

    return dat_davg, dat_ts, time_nm

def do_fft(dat, n):
    # carries out fast fourier transform method of checking for periocity in time series
    fft_result = fft(np.asarray(dat))
    power_spec = np.abs(fft_result) ** 2
    freq = np.fft.fftfreq(n, 1) # Frequency axis (in cycles per data point ie time step)
    return freq, power_spec
def get_meshmask(path_mm, file_mm,row_j=175,col_i=75):
    # returns masking arrays from the meshmask given a row,col input
    # used for depth integration

    with nc.Dataset(os.path.join(path_mm, file_mm), 'r') as mesh:
        tmask = mesh.variables['tmask'][:]
        #e3t0 the depth bin widths
        e3t0 = mesh.variables['e3t_0'][:]
        # the centres of the depth bins
        gdept_0 = mesh.variables['gdept_0'][:]

    gdept_0 = gdept_0[0, :, row_j, col_i]
    e3t0 = e3t0[0, :, row_j, col_i]
    tmask = tmask[0, :, row_j, col_i]

    return(tmask,gdept_0,e3t0)

def depth_int(data_in, depth_min, depth_max, gdept_0, e3t0):
    # calculates depth integrated avg account for variable depth widths
    # Note this assumes that data and tmask have same depth strata
    # input -
    #  data_in - a numpy array or xr DataArray, must contain deptht
    #  depth_min - the minimum depth
    #  depth_max - the max depth
    #  gdept_0 - from NEMO model mesh mask, centre of vert depth bins
    #  e3t0 - the widths of depth bins (varies)
    #
    # note could be improved b/c avg being taken from closest g_dep point
    # even if that is smaller than dep min or larger than dep max
    # instead we could re-interpolate to the desired dep min max first

    # find closest idx in depth to the min / max depths
    diffs_min = np.abs(gdept_0 - depth_min)
    min_idx = np.argmin(diffs_min)
    diffs_max = np.abs(gdept_0 - depth_max)
    max_idx = np.argmin(diffs_max)
    e3t0_2 = e3t0[min_idx:max_idx + 1]

    # !this assumes data and mask have same depth strata!
    data_in_trim = data_in[min_idx:max_idx + 1]

    # replace values in e3t where obs are nans with nans
    nan_indices = np.isnan(data_in_trim)

    # add dimensions to e3t to match shape of data
    dat_shape = data_in_trim.shape
    e3t0_3 = e3t0_2[:, np.newaxis]
    e3t0_3 = np.repeat(e3t0_3, dat_shape[1], axis=1)

    e3t0_3[nan_indices] = np.nan

    # average using depth widths as weights
    sum_e3t = np.nansum(e3t0_3, axis=0)
    e3t0_weights = np.where(sum_e3t != 0, e3t0_3 / sum_e3t[np.newaxis, :], 0)
    #crosscheck_weights = np.nansum(e3t0_weights, axis=0)
    #print(crosscheck_weights) # should sum to one or be zeros (zeros replaced with nan later)
    weighted_vals = data_in_trim * e3t0_weights
    dat_avg = np.sum(weighted_vals, axis=0)  # sum instead of nansum keeps it as xarray dataarray
    dat_avg = dat_avg.where(dat_avg != 0, np.nan)

    return(dat_avg)