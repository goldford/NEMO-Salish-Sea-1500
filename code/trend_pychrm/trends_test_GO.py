import os
import numpy as np
import xarray as xr
from GO_tools import get_meshmask, depth_int, do_fft, do_seasonal_avg
from scipy.stats import linregress
from linregress_GO import linregress_GO
from scipy import stats
import datetime as dt
import mannkendall as mk
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf
from statsmodels.graphics.tsaplots import plot_acf
from scipy.fft import fft, ifft

# ==== Paths ====
meshm_p = '../../data/mesh mask/'
meshm_f = 'mesh_mask_20210406.nc'
#dat_f = 'Nanoose_obs_anom_temp_1970-2005.nc'
dat_f = 'Nanoose_obs_anom_temp_1980-2018.nc'
#dat_f = 'Nanoose_obs_anom_temp_1970-2018.nc'
#dat_f = 'RUN216mod_anom_temp_1980-2018.nc'
#dat_f = 'RUN203mod_anom_temp_1980-2018.nc'

# ==== Spatial-temporal settings ====
depth_min = 0
depth_max = 30

# ==== Load Data ====
var = "temperature"
dat = xr.open_dataset(os.path.join('../climatol_intermediate_files/',dat_f))
dat = dat[var]

# ==== Depth integration ====
tmask, gdept_0, e3t0 = get_meshmask(meshm_p,meshm_f)
dat_davg = depth_int(dat, depth_min, depth_max, gdept_0, e3t0)

# ==== Deal with missing values ====
remove_nans = True
if remove_nans:
    #nan_mask = np.logical_or(np.isnan(dat_davg), np.isnan(time_numeric))
    nan_mask = np.isnan(dat_davg)
    dat_davg = dat_davg[~nan_mask]

# ==== Time Averaging ====
# monthly, seasonal, and annual come in handy later
# assumes dataset is biweekly
time_inc = "seasonal" # biweekly, monthly, seasonal, annual
print("Using ", time_inc, " avg")
dat_davg, dat_ts, time_nm = do_seasonal_avg(dat_davg,time_inc)

# ==== Deal with time, datetime, etc ====
start_date = dat_davg[time_nm].values[0]
#  assumes 'time_counter' is in datetime64 format
time_numeric = (dat_davg[time_nm] - start_date) / np.timedelta64(1, 'D')
time_datetime64 = np.array(dat_davg[time_nm], dtype='datetime64[s]')
#  painful to convert to py datetime but is required for mk
time_pydatetime = np.array([dt.datetime(year, month, day) for year, month, day in zip(dat_davg[time_nm].dt.year.values,
                                                                                dat_davg[time_nm].dt.month.values,
                                                                            dat_davg[time_nm].dt.day.values)])

# ==== Estimate trend ====
LR_tf = True # Linear Regression
if LR_tf:
    slope_LR, inter_LR, rval_LR, pval_LR, _ = linregress(time_numeric.values, dat_davg.values)
    slope_LR_yr = slope_LR * 365 # returns slope in /day (the units of time_numeric
    inter_LR_yr = inter_LR * 365  # returns slope in /day (the units of time_numeric
    print("Slope from LR: ", slope_LR_yr)
    print("p_val from LR (no adjustment for ESS): ", pval_LR)

SS_tf = True # Sen's slope w/ CI
resolution = 0.0001# sig digits for ties sensitivity
if SS_tf:
    dat_davg_np = np.asarray(dat_davg)
    t = mk.mkt.nb_tie(dat_davg_np, resolution)
    (s, n) = mk.mks.s_test(dat_davg_np, time_pydatetime)
    k_var = mk.mkt.kendall_var(dat_davg_np, t, n)
    z = mk.mks.std_normal_var(s, k_var)
    mk_out = mk.mks.sen_slope(time_pydatetime, dat_davg_np,k_var) # units of second!
    mult = 3600 * 24 * 365.25
    slope_yr = mk_out[0] * mult; slope_yr_min = mk_out[1] * mult; slope_yr_max = mk_out[2] * mult
    print("Sens slope, no seas: ", slope_yr , " Min: ", slope_yr_min, " Max: ", slope_yr_max)


# ==== Detrend and test normality of resids ====
# (to do) - try detrend seasonally, use sen's slope
dat_davg_dtrnd = dat_davg - (slope_LR * time_numeric + inter_LR)

# Test resid for norm using Shapiro-Wilks
residuals = dat_davg_dtrnd
stat, p_value = stats.shapiro(residuals)
print(f'Statistic: {stat}, p-value: {p_value}')

if p_value > 0.05:
    print("Residuals likely are normally distributed. Parametric methods could be used for CI estimation.")
    type_CI = "parametric"
else:
    print("Residuals are unlikely to be drawn from a normal distribution. Non-parametric methods for CI estimation should be considered.")
    type_CI = "non-parametric"

# To do - are they non-normal if averaged to monthly, seasonal, annual?

# ==== test for seasonality using FFT ====
freq, power_spec = do_fft(dat_davg_dtrnd, len(dat_davg_dtrnd))
peak_indices = np.argsort(power_spec)[::-1][:5]  # Top 5 peaks
peak_frequencies = np.abs(freq[peak_indices])    # corresponding frequencies
peak_freq_period_yrs = 1/peak_frequencies*1/dat_ts

# plot power spec
len_ts = len(dat_davg_dtrnd)
plt.plot(freq[:len_ts // 2], power_spec[:len_ts // 2]) # Exclude the negative freqs
plt.title('Power Spectrum')
plt.xlabel('Frequency (cycles per data point)')
plt.ylabel('Power')
plt.show()

# check for sub-annual seasonal periodicities
print("Top peak frequencies in years:", peak_freq_period_yrs)
print(np.any(peak_freq_period_yrs <= 1))
subannual_freq_tf = np.any(peak_freq_period_yrs <= 1)
if subannual_freq_tf:
    print("Subannual cycles are detected in anomaly residuals.")
else:
    print("Subannual cycles are not detected in anomaly residuals.")

type_CI = "parametric"
if type_CI == "parametric":

    print("proceeding with CI estimation using decorrelation time scale")
    # - ACF and Decorr Time Scale (AR1)
    alpha = 0.05
    nlags = 48

    acf_values = acf(residuals.values, adjusted=True, fft=False, alpha=alpha, nlags=nlags, missing="conservative")
    plot_acf(residuals.data, adjusted=True, fft=False, alpha=alpha, lags=len(residuals)-1)
    plt.title('Autocorrelation Function with Confidence Intervals')
    plt.xlabel('Lag')
    plt.ylabel('Autocorrelation (-1, 1)')
    plt.show()

    AR1_coef = acf_values[0][1]
    DTS = (1 + AR1_coef) / (1 - AR1_coef) # decorr. time scale in ts units
    ESS = len(dat_davg) / DTS          # effective sample size

    print("The AR(1) Coefficient:", AR1_coef)
    print("The Decorrelation Time Scale (DTS; units: yr): ", DTS / dat_ts)
    print("The Effective Sample Size (ESS): ", ESS)

    # adjust sample size downward
    slope_dst, inter_dst, r_val_dst, p_val_dst, std_err_dst, _ = linregress_GO(time_numeric.values,
                                                                  dat_davg.values, ESS)
    df = ESS - 2
    t_critical = stats.t.ppf((1 + (1-alpha)) / 2, df)
    margin_of_error = t_critical * std_err_dst
    confidence_interval = ((slope_dst - margin_of_error) * 365, (slope_dst + margin_of_error) * 365)
    print(f"Slope: {slope_dst*365}")
    print(f"Confidence Interval (95%): {confidence_interval}")
    print(f"P val:  {p_val_dst}")

# 2) Non-Parametric Trend Detection
# - MK w/ Prewhitening (e.g., 3PW)
type_CI = "nonparametric"
resolution = 0.0001 # of measurements to use in ties calc in Kendalls Tau and S
if time_inc == "annual": seasons_per_year = 1
elif time_inc == "seasonal": seasons_per_year = 4
elif time_inc == "monthly": seasons_per_year = 12
elif time_inc == "biweekly": seasons_per_year = 24

if type_CI == "nonparametric":
    print("MK w/ Sens and 3PW method")
    seasonal_dates = []
    dat_multi = []
    if seasons_per_year == 1:
        out = mk.mk_temp_aggr(time_pydatetime, np.asarray(dat_davg), resolution)
        print(out)
    else:
        for season in range(1, seasons_per_year + 1):
            if seasons_per_year == 24:
                indices = np.where(
                    (dat_davg[time_nm].dt.dayofyear >= (season * 15) - 10) &
                    (dat_davg[time_nm].dt.dayofyear < ((season + 1) * 15) - 10))[0]
            elif seasons_per_year == 12:
                indices = np.where(dat_davg[time_nm].dt.month == season)
            elif seasons_per_year == 4:
                if season == 1: season_name = "winter"
                elif season == 2: season_name = "spring"
                elif season == 3: season_name = "summer"
                elif season == 4: season_name = "fall"
                else: print("error w/ seasons")
                indices = np.where(dat_davg['season'] == season_name)
            seasonal_dates.append((time_pydatetime[indices]))
            dat_multi.append((dat_davg[indices]))
        out = mk.mk_temp_aggr(seasonal_dates, dat_multi, resolution)
        for n in range(seasons_per_year):
            print('Season {ind}:'.format(ind=n + 1), out[n])
        print('Overall: ', out[len(out)-1])


# to do, bootstrap

# else - to do