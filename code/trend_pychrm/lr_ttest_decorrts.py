# Created Jan 19 by G Oldford
# Purpose: conduct a ordinary least squared regression on the Nanoose Dataset
#          using a t-test that adjusts for autocorrelation using the decorrelation time scale
#
# Context: the mannkendall and pymannkendall libraries both provided challenges
#          the most important of which was apparent bug in mannkendall
#          and lack of clarity about seasonal analysis results https://github.com/mannkendall/Python/issues

import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import scipy
from scipy.stats import linregress, ttest_1samp
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
import os
import statsmodels
from statsmodels.tsa.stattools import acf
from statsmodels.graphics.tsaplots import plot_acf
print(statsmodels.__version__)
import sys
print(sys.version)
from linregress_GO import linregress_GO
from GO_tools import get_meshmask, depth_int


meshm_p = '../../data/mesh mask/'
meshm_f = 'mesh_mask_20210406.nc'
tmask, gdept_0, e3t0 = get_meshmask(meshm_p,meshm_f)

# Open the netCDF file
# Replace 'your_file.nc' with the actual file path
#obs_f = 'Nanoose_obs_anom_temp_1970-2005.nc'
obs_f = 'Nanoose_obs_anom_temp_1980-2018.nc'
#obs_f = 'Nanoose_obs_anom_temp_1970-2018.nc'
obs = xr.open_dataset(os.path.join('../climatol_intermediate_files/',obs_f))
mod_f = 'RUN216mod_anom_temp_1980-2018.nc'
#mod_f = 'RUN203mod_anom_temp_1980-2018.nc'
mod = xr.open_dataset(os.path.join('../climatol_intermediate_files/',mod_f))

time_counter_obs = obs['time_counter']
time_datetime_obs = np.array(time_counter_obs, dtype='datetime64[s]')
time_counter_mod = mod['time_counter']
time_datetime_mod = np.array(time_counter_mod, dtype='datetime64[s]')

# data's time increment in years
dat_ts = 24

var = 'temperature'
var_dat_obs = obs[var]
var_dat_mod = mod[var]
depth_min = 4.5
depth_max = 400
# get depth strata widths and depth integrate
var_dat_avg_obs = depth_int(var_dat_obs, depth_min, depth_max, gdept_0, e3t0)
var_dat_avg_mod = depth_int(var_dat_mod, depth_min, depth_max, gdept_0, e3t0)

start_date_obs = time_counter_obs.values[0]
start_date_mod = time_counter_mod.values[0]

# Assuming 'time_counter' is in datetime64 format, convert it to numeric values
# representing time in days since a reference point (e.g., 1970-01-01)
time_numeric_obs = (time_counter_obs - start_date_obs) / np.timedelta64(1, 'D')
time_numeric_mod = (time_counter_mod - start_date_mod) / np.timedelta64(1, 'D')

# linear regression
# Check for NaN values in variable and time_numeric
nan_mask_obs = np.logical_or(np.isnan(var_dat_avg_obs), np.isnan(time_numeric_obs))
nan_mask_mod = np.logical_or(np.isnan(var_dat_avg_mod), np.isnan(time_numeric_mod))

var_dat_avg_obs_nonan = var_dat_avg_obs[~nan_mask_obs]
var_dat_avg_mod_nonan = var_dat_avg_mod[~nan_mask_mod]
time_numeric_obs_nonan = time_numeric_obs[~nan_mask_obs]
time_numeric_mod_nonan = time_numeric_mod[~nan_mask_mod]
time_datetime_obs_nonan = time_datetime_obs[~nan_mask_obs]
time_datetime_mod_nonan = time_datetime_mod[~nan_mask_mod]

slope_o, inter_o, r_val_o, p_val_o, std_err_o = linregress(time_numeric_obs_nonan.values, var_dat_avg_obs_nonan.values)
slope_m, inter_m, r_val_m, p_val_m, std_err_m = linregress(time_numeric_mod_nonan.values, var_dat_avg_mod_nonan.values)

# Print regression results

print("Obs Slope:", slope_o)
print("Mod Slope:", slope_m)
print("Obs Intercept:", inter_o)
print("Mod Intercept:", inter_m)
print("Obs R-squared:", r_val_o**2)
print("Mod R-squared:", r_val_m**2)
print("Obs P-value:", p_val_o)
print("Mod P-value:", p_val_m)
print("Obs Standard Error:", std_err_o)
print("Mod Standard Error:", std_err_m)

# decorrelation time scale
# Assuming 'data' is your time series data
alpha = 0.05

# choice here is to use or not use nans in acf
notrend_obs = var_dat_avg_obs_nonan - (slope_o*time_numeric_obs_nonan.values+inter_o)
notrend_mod = var_dat_avg_mod_nonan - (slope_m*time_numeric_mod_nonan.values+inter_m)
# notrend_obs = var_dat_avg_obs - (slope_o*var_dat_avg_obs+inter_o)
# notrend_mod = var_dat_avg_mod - (slope_m*var_dat_avg_mod+inter_m)

# my understanding is adjusted = T used if mean has been removed, which is climatol process
# see Thompson et al 2014 Chp 5, p. 429
acf_values_obs = acf(notrend_obs.values, adjusted=True, fft=False, alpha=alpha, nlags=96, missing="conservative")
acf_values_mod = acf(notrend_mod.values, adjusted=True,fft=False, alpha=alpha, nlags=96, missing="conservative")
#print("Obs ACF Values:", acf_values_obs)
#print("Mod ACF Values: ", acf_values_mod)

# since lag-0 acf would equate to 1, the first index here is lag-1 (idx=1)
AR1_coef_obs = acf_values_obs[0][1];
AR1_coef_mod = acf_values_mod[0][1]
print("The AR(1) Coefficient, Obs:", AR1_coef_obs)
print("The AR(1) Coefficient, Mod:", AR1_coef_mod)
# the decorrelation time scale may be calculated a number of diff ways.
# See Storch & Zwiers, 1999, p. 371
# For AR(1) process, the DTS can be found simply, using lag-1 AR coeff
DTS_obs = (1 + AR1_coef_obs) / (1 - AR1_coef_obs) * 1/dat_ts
DTS_mod = (1 + AR1_coef_mod) / (1 - AR1_coef_mod) * 1/dat_ts
print("The DTS (yr), Obs: ", DTS_obs)
print("The DTS (yr), Mod: ", DTS_mod)
# sanity check:
# DTS_obs_2 = 1 + 2 * (1 + AR1_coef_obs)
# print("Sanity check - The DTS, Obs: ", DTS_obs)
# effective sample size (ESS)
ESS_obs = (len(var_dat_avg_obs_nonan)*1/dat_ts) / DTS_obs
ESS_mod = (len(var_dat_avg_mod_nonan)*1/dat_ts) / DTS_mod

# now we adjust sample size downward

slope_o, inter_o, r_val_o, p_val_o, std_err_o, _ = linregress_GO(time_numeric_obs_nonan.values,
                                                              var_dat_avg_obs_nonan.values,
                                                              ESS_obs)
slope_m, inter_m, r_val_m, p_val_m, std_err_m, _ = linregress_GO(time_numeric_mod_nonan.values,
                                                              var_dat_avg_mod_nonan.values,
                                                              ESS_obs)
print("Obs P-value:", p_val_o)
print("Mod P-value:", p_val_m)

# Plot ACF with confidence interval bounds
plot_acf(notrend_obs.data, adjusted=True, fft=False, alpha=alpha, lags=96)
plt.title('Obs. Autocorrelation Function with Confidence Intervals')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.show()

plot_acf(notrend_mod.data, adjusted=True, fft=False, alpha=alpha, lags=120)
plt.title('Model Autocorrelation Function with Confidence Intervals')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.show()

# pre-calculated and hard-coded based on interpretation of plots above

# Plot the data and regression line
fig, ax = plt.subplots()

# Plot the cumulative difference
nan_mask = np.logical_or(np.logical_or(np.isnan(var_dat_avg_obs), np.isnan(time_numeric_obs)),
                         np.logical_or(np.isnan(var_dat_avg_mod), np.isnan(time_numeric_mod))
                         )
var_dat_avg_obs_nonan2 = var_dat_avg_obs[~nan_mask]
var_dat_avg_mod_nonan2 = var_dat_avg_mod[~nan_mask]
time_numeric_obs_nonan2 = time_numeric_obs[~nan_mask]
time_numeric_mod_nonan2 = time_numeric_mod[~nan_mask]
time_datetime_obs_nonan2 = time_datetime_obs[~nan_mask]
time_datetime_mod_nonan2 = time_datetime_mod[~nan_mask]

#print([0] + [sum(var_dat_avg_obs[:i + 1]) - sum(var_dat_avg_mod[:i + 1]) for i in range(len(time_datetime_obs) - 1)])
ax.fill_between(time_datetime_obs_nonan2, var_dat_avg_obs_nonan2, var_dat_avg_mod_nonan2, color='grey', label='', alpha=0.7, zorder=0)

ax.scatter(time_datetime_obs, var_dat_avg_obs, label='Obs', s=1, color='r', marker='o')
ax.scatter(time_datetime_mod, var_dat_avg_mod, label='Mod', s=1, color='g', marker='x')

ax.plot(time_datetime_obs, inter_o + slope_o * mdates.date2num(time_datetime_obs), 'r', label='')
ax.plot(time_datetime_mod, inter_m + slope_m * mdates.date2num(time_datetime_mod), 'g', label='')

# Format the x-axis as dates
ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

# Annotate the plot with the slope value
slope_label = "Obs Slope: {:.4f} /yr{}".format(slope_o * 365, "\n") + "Mod Slope: {:.4f} /yr".format(slope_m * 365)

ax.text(0.5, 0.05, slope_label, transform=ax.transAxes, fontsize=10, color='black')

plt.title("Run WITH Bias Correction, Depths " + str(depth_min) + " to " + str(depth_max) + " m")

plt.xlabel('Time')
plt.ylabel('Temperature Anomaly')
plt.legend()
plt.show()

