
# Created Jan 2024 by G Oldford
# purpose: generate trend statistics for model and obs data from the Nanoose Station time series
# input:
#    meshm_p     - string; path to mesh mask from NEMO which contains depths and coords
#    meshm_f     - string; mesh mask file name
#    dat_f       - string; netcdf data file of observations, assumed biweekly and spatially averaged (dims: depth, time)
#    depth_min   - float; minimum depth to use (during depth integration)
#    depth_max   - float; maximum depth to use (during depth integration)
#    var         - string; the variable being analysed (e.g., "temperature") - only tested with temp for now
#    resolution  - float; the precision of the measurement (sig digits) for ties calculated during some stats calcs
#    remove_nans - True/False;
#    time_inc    - string; sub-annual time increment for temporal avg ("biweekly", "monthly", "seasonal", "annual")
#    alpha_CI    - the alpha confidence interval threshold (e.g., 0.05 for 95%)
# returns:
#    slopes
#    summary_stats_all
#    DTS_CI_out        - results of trend test using decorr. time scale
#                      - list of dict elements; {season: {slope, inter, r_val, p_val, std_err, slope_lcl, slope_ucl}}
#    mk_3pw_out        - results of '3PW' trend test using Mann-Kendall, see https://doi.org/10.5194/amt-13-6945-2020
#                      - list of dict elements; {P= max(P_PW, P_TFPW_Y);'ss' (float): statistical significance:
#                          alpha_MK if the test is ss at the alpha confidence level. Defaults to 95.
#                          0 if the test is not ss at the alpha_MK confidence level.
#                          -1 if the test is a TFPW_Y false positive at alpha_MK confidence level
#                          -2 if the test is a PW false positive at alpha_MK confidence level
#                          'slope' (float): Sen's slope in units/y
#                           'ucl' (float): upper confidence level in units/y
#                           'lcl' (float): lower confidence level in units/y
# notes:
#    - the modelled data inputs may or may not be 'blind' (nans where obs has gaps)
#    - for time averaging, the finest resolution assumed present in the data is biweekly
#    - to do: bootstrap, make below a function that can be looped, auto-select trend test after stats done,
#             introduce some kind of harmonic analysis or ESD


import os
import numpy as np
import xarray as xr
from GO_tools import (get_meshmask, depth_int, do_fft, do_seasonal_avg,
                      do_ss_seasons, do_lr_seasons, prep_data_seas, do_summary_stats,
                      dst_trend_sig, mk3pw_trend_sig)
import datetime as dt
from tabulate import tabulate

# ==== Paths ====
meshm_p = '../../data/mesh mask/'
meshm_f = 'mesh_mask_20210406.nc'
dat_p = '../climatol_intermediate_files/'
#dat_f = 'Nanoose_obs_anom_temp_1970-2005.nc'
dat_f = 'Nanoose_obs_anom_temp_1980-2018.nc'
#dat_f = 'Nanoose_obs_anom_temp_1970-2018.nc'
#dat_f = 'RUN216mod_anom_temp_1980-2018.nc'
#dat_f = 'RUN203mod_anom_temp_1980-2018.nc'

# ==== Params for Analysis ====
depth_min = 4.5
depth_max = 400
var = 'temperature'
resolution = 0.0001 # measurement precision, for Kendall's tau etc
remove_nans = True
alpha_DTS_CI = 0.05
alpha_MK = 0.95
time_inc = 'biweekly' # biweekly, monthly, seasonal, annual
resolution = 0.0001  # of measurements to use in ties calc in Kendalls Tau and S
if time_inc == 'annual': seasons_per_year = 1
elif time_inc == 'seasonal': seasons_per_year = 4
elif time_inc == 'monthly': seasons_per_year = 12
elif time_inc == 'biweekly': seasons_per_year = 24

# ==== Load Data ====
dat = xr.open_dataset(os.path.join(dat_p,dat_f))
dat = dat[var]

# ==== Depth integration ====
tmask, gdept_0, e3t0 = get_meshmask(meshm_p,meshm_f)
dat_davg = depth_int(dat, depth_min, depth_max, gdept_0, e3t0)

# ==== Deal with missing values ====
if remove_nans:
    #nan_mask = np.logical_or(np.isnan(dat_davg), np.isnan(time_numeric))
    nan_mask = np.isnan(dat_davg)
    dat_davg = dat_davg[~nan_mask]

# ==== Time Averaging ====
print("Using ", time_inc, " avg")
dat_davg, dat_ts, time_dim = do_seasonal_avg(dat_davg,time_inc)

# ==== Deal with time, datetime, etc ====
start_date = dat_davg[time_dim].values[0]
#  assumes 'time_counter' is in datetime64 format
dat_numerictime = (dat_davg[time_dim] - start_date) / np.timedelta64(1, 'D')
time_datetime64 = np.array(dat_davg[time_dim], dtype='datetime64[s]')
#  painful to convert to py datetime but is required for mk
dat_pytime = np.array([dt.datetime(year, month, day) for year, month, day in zip(dat_davg[time_dim].dt.year.values,
                                                                                dat_davg[time_dim].dt.month.values,
                                                                            dat_davg[time_dim].dt.day.values)])
# split seasons
dat_seas, dat_pytime_seas, dat_timenumeric_seas = prep_data_seas(dat_davg, dat_pytime, dat_numerictime,
                                                                  seasons_per_year, time_dim)
# slopes from linear regression and alt method theil-sen aka sen's
LR_slopes = do_lr_seasons(dat_davg, dat_seas, dat_numerictime, dat_timenumeric_seas, seasons_per_year)
SS_slopes = do_ss_seasons(dat_davg, dat_seas, dat_pytime, dat_pytime_seas, seasons_per_year, resolution)

# ==== Detrend ====
# as in mannkendall lib (mk_white.py)
slope_method = 'SS' # To do - make "LR" work
if slope_method == 'SS': slopes = SS_slopes
elif slope_method == "LR": slopes = LR_slopes
else: print("error with slope method")
# note that low p-vals in Het-White and Goldfeld-Quandt tests suggests evidence against homoscedasticity
summary_stats_all = do_summary_stats(slopes, slope_method, dat_ts, dat_davg, dat_numerictime, dat_seas, dat_timenumeric_seas)

# ====  Parametric Trend Detection =====
# use LR or SS with adjusted ESS from DST
# to-do: make this automatic choice
type_CI = 'parametric'
if type_CI == 'parametric':
    DTS_CI_out = dst_trend_sig(summary_stats_all, dat_ts, alpha_DTS_CI, dat_davg, dat_seas, dat_numerictime, dat_timenumeric_seas)

# ==== Non-Parametric Trend Detection ====
# - MK w/ Prewhitening (e.g., 3PW)
type_CI = 'nonparametric'
if type_CI == 'nonparametric':
    print("MK w/ Sens and 3PW method")
    mk_3pw_out = mk3pw_trend_sig(seasons_per_year, resolution, alpha_MK, dat_davg, dat_pytime, dat_seas, dat_pytime_seas)
print(mk_3pw_out)

# write to tables
headers = list(summary_stats_all[0]['season 1'].keys())
rows = []
seas = 1
for s in summary_stats_all:
    if seas == len(summary_stats_all):
        seas_name = 'all seasons'
        row = s[seas_name].values()
    else:
        seas_name = 'season ' + str(seas)
        row = s[seas_name].values()
    seas += 1
    row = [seas_name] + list(row)
    rows.append(row)


table = tabulate(rows, headers, tablefmt='grid')
print(table)
# to do

# temp, to get a plot
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
dat_f = 'RUN216mod_anom_temp_1980-2018.nc'
dat = xr.open_dataset(os.path.join(dat_p,dat_f))
dat = dat[var]
mod_davg = depth_int(dat, depth_min, depth_max, gdept_0, e3t0)
if remove_nans:
    #nan_mask = np.logical_or(np.isnan(dat_davg), np.isnan(time_numeric))
    nan_mask = np.isnan(mod_davg)
    mod_davg = mod_davg[~nan_mask]
ax.fill_between(time_datetime_obs_nonan2, dat_davg, mod_davg, color='grey', label='', alpha=0.7, zorder=0)
plt.title("Run WITH Bias Correction, Depths " + str(depth_min) + " to " + str(depth_max) + " m")
plt.xlabel('Time')
plt.ylabel('Temperature Anomaly')
plt.legend()
plt.show()
