# created Jan 2026 by G Oldford
# purpose: analyze trends at Nanoose, model and obs

# to do:
#    test if results different if not using anomalies

import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import datetime as dt
import mannkendall as mk
import netCDF4 as nc
import warnings




# annoying but necessary conversion from datetime64[ns] to py dt.dateime
def dt64_to_pydt(data_in,tvar,month=1):
    # created by G Oldford Jan 2024
    # assumes a tvar in xr is np datetime64
    # and converts to python datetime
    # input
    #   data_in - xarray DataArray
    #   tvar - string ref to xarray DataArray time variable
    # returns

    # if averaging annually (date is year only), we still need day
    # and month for pd.datetime required by mannkendall lib
    if tvar == 'year':
        years = np.asarray(data_in[tvar])
        time_np_dt = np.array([dt.datetime(year, month, 1) for year in years])
    else:
        years = data_in[tvar].dt.year
        months = data_in[tvar].dt.month
        days = data_in[tvar].dt.day
        time_np_dt = np.array([dt.datetime(year, month, day) for year, month, day in zip(years.values, months.values, days.values)])

    return(time_np_dt)

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

meshm_p = '../../data/mesh mask/'
meshm_f = 'mesh_mask_20210406.nc'
tmask,gdept_0,e3t0 = get_meshmask(meshm_p,meshm_f)

#

variables_dict = {'temperature'}

period_dict = {"MC07: 1970 to 2005":{"year_min":1970, "year_max": 2005,
                                     "file_anom_temp": "Nanoose_obs_anom_temp_1970-2005.nc",
                                     "file_raw_temp": "obs_temperature_1970-2005-bimonthly_timeseries_mc07.nc"},
                   "HOTSS period: 1980 to 2018":{"year_min":1980, "year_max":2018,
                                          "file_anom_temp": "Nanoose_obs_anom_temp_1980-2018.nc",
                                          "file_raw_temp": "obs_temperature_1980-2018-bimonthly_timeseries.nc"},
                   "Nanoose Full: 1970 to 2018":{"year_min":1970, "year_max":2018,
                                          "file_anom_temp": "Nanoose_obs_anom_temp_1970-2018.nc",
                                          "file_raw_temp": "obs_temperature_1970-2018-bimonthly_timeseries.nc"},
                   }

depth_lev_dict = {"4.5to400":{"depth_min":4.5, "depth_max": 400},
                   "0to30":{"depth_min":0, "depth_max":30},
                   "30to150":{"depth_min":30, "depth_max":150},
                   "150anddeeper":{"depth_min":150, "depth_max":400}
                   }
seasons = {
    'winter': [12, 1, 2],
    'spring': [3, 4, 5],
    'summer': [6, 7, 8],
    'fall': [9, 10, 11]
}

obs_f = 'Nanoose_obs_anom_temp_1970-2005.nc'
#obs_mc07_f = 'obs_temperature_1970-2005-bimonthly_timeseries_mc07.nc'
obs = xr.open_dataset(os.path.join('../climatol_intermediate_files/',obs_f))

variable = 'temperature'
for d_lev in depth_lev_dict.keys():
    print("Depth level: " + d_lev)
    depth_min = depth_lev_dict[d_lev]['depth_min']
    depth_max = depth_lev_dict[d_lev]['depth_max']

    # get depth integrated average
    obs_avg = depth_int(obs[variable], depth_min, depth_max, gdept_0, e3t0)
    obs_avg_dt_py = dt64_to_pydt(obs_avg['time_counter'], 'time_counter')

    # THESE ARE WRONG - NEED PROPER DATES can't have seasons with same dates

    # all-together seasonal analysis using the seas. procedure in mannkendall
    # formats data as required by mannkendall, treating the time step in prepped data as 'season'
    seasons_per_year = 24
    print("Method 1a: seasonal analysis using " + str(seasons_per_year) + " mannkendall procedure")
    obs_avg_dt_seas = []
    obs_avg_multi = []
    for seas in range(1, seasons_per_year + 1):
        indices = np.where(
            (obs_avg['time_counter'].dt.dayofyear >= (seas * 15) - 10) &
            (obs_avg['time_counter'].dt.dayofyear < ((seas + 1) * 15) - 10))[0]
        obs_avg_dt_seas.append(np.asarray(obs_avg_dt_py[indices]))
        obs_avg_multi.append(np.asarray(obs_avg[indices]))

    with warnings.catch_warnings():
        #warnings.simplefilter("ignore", category=UserWarning)
        out = mk.mk_temp_aggr(obs_avg_dt_seas, obs_avg_multi, 0.00001)
        for n in range(seasons_per_year):
            print('Season {ind}:'.format(ind=n + 1), out[n])
        print('Combined yearly trend: ', out[seasons_per_year])

    seasons_per_year = 12
    print("Method 1b: seasonal analysis using " + str(seasons_per_year) + " mannkendall procedure")
    obs_avg_dt_seas = []
    obs_avg_multi = []
    for seas in range(1, seasons_per_year + 1):
        indices = np.where(
            obs_avg['time_counter'].dt.month == seas)
        month = seas
        obs_avg_seas = obs_avg[indices].groupby("time_counter.year").mean(dim="time_counter")
        obs_avg_seas_dt_py = dt64_to_pydt(obs_avg_seas['year'], 'year', month)
        obs_avg_multi.append(np.asarray(obs_avg_seas))
        obs_avg_dt_seas.append(np.asarray(obs_avg_seas_dt_py))

    with warnings.catch_warnings():
        #warnings.simplefilter("ignore", category=UserWarning)
        out = mk.mk_temp_aggr(obs_avg_dt_seas, obs_avg_multi, 0.00001)
        for n in range(seasons_per_year):
            print('Season {ind}:'.format(ind=n + 1), out[n])
        print('Combined yearly trend: ', out[seasons_per_year])


    print("Method 1c: increase granularity to seasonal before analysis using mk seasonal")
    obs_avg_dt_py_multi = []
    obs_avg_multi = []
    for season_name, season_months in seasons.items():
        indices = np.where(obs_avg['time_counter'].dt.month.isin(season_months))
        obs_avg_seas = obs_avg[indices].groupby("time_counter.year").mean(dim="time_counter")
        obs_avg_multi.append(np.asarray(obs_avg_seas))
        obs_avg_seas_dt_py = dt64_to_pydt(obs_avg_seas['year'], 'year', season_months[1])
        obs_avg_dt_py_multi.append(obs_avg_seas_dt_py)

    with warnings.catch_warnings():
        #warnings.simplefilter("ignore", category=UserWarning)
        out = mk.mk_temp_aggr(obs_avg_dt_py_multi, obs_avg_multi, 0.00001)
        for n in range(len(seasons)):
            print('Season {ind}:'.format(ind=n + 1), out[n])
        print('Combined yearly trend: ', out[len(seasons)])

    print("Sen's Slope w/ MannKendall and 3PW Procedure:")
    # annual average analysis
    print('Method 2: annually averaging within season and truncating to season before mk test')
    for season_name, season_months in seasons.items():
        print('Season: ' + season_name)

        season_data = obs_avg.sel(time_counter=obs_avg.time_counter.dt.month.isin(season_months))
        season_data_ann = season_data.groupby("time_counter.year").mean(dim="time_counter")
        seasonal_data_ann_np = np.asarray(season_data_ann)
        season_data_ann_dt_py = dt64_to_pydt(season_data_ann['year'], 'year', season_months[1])

        # Catch and suppress UserWarnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            out = mk.mk_temp_aggr(season_data_ann_dt_py, seasonal_data_ann_np, 0.00001)

        # not sure why two identical elements returned, so just print one
        if len(out)>1:
            for n in range(len(out)-1):
                print(out[n])
        else:
            print(out)


# to do: take seasonal average, analyze
# method for decorrelation time scale
