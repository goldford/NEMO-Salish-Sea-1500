import os
import numpy as np
import xarray as xr
from GO_tools import (get_meshmask, depth_int, do_fft, do_seasonal_avg,
                      do_ss_seasons, do_lr_seasons, prep_data_seas, do_summary_stats,
                      dst_trend_sig, mk3pw_trend_sig)
import datetime as dt

# below to be removed
from linregress_GO import linregress_GO
from scipy import stats
import mannkendall as mk
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.stats.diagnostic import het_goldfeldquandt, het_white
from statsmodels.tools.tools import add_constant

# ==== Paths ====
meshm_p = '../../data/mesh mask/'
meshm_f = 'mesh_mask_20210406.nc'
#dat_f = 'Nanoose_obs_anom_temp_1970-2005.nc'
dat_f = 'Nanoose_obs_anom_temp_1980-2018.nc'
#dat_f = 'Nanoose_obs_anom_temp_1970-2018.nc'
#dat_f = 'RUN216mod_anom_temp_1980-2018.nc'
#dat_f = 'RUN203mod_anom_temp_1980-2018.nc'

# ==== Params for Analysis ====
depth_min = 4.5
depth_max = 400
var = "temperature"
resolution = 0.0001 # measurement precision, for Kendall's tau etc
remove_nans = True
time_inc = "seasonal" # biweekly, monthly, seasonal, annual
if time_inc == "annual": seasons_per_year = 1
elif time_inc == "seasonal": seasons_per_year = 4
elif time_inc == "monthly": seasons_per_year = 12
elif time_inc == "biweekly": seasons_per_year = 24

# ==== Load Data ====
dat = xr.open_dataset(os.path.join('../climatol_intermediate_files/',dat_f))
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
# monthly, seasonal, and annual come in handy later (assumes dataset is biweekly)
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
# (to do) - try detrend seasonally
# as in mannkendall lib (mk_white.py), use sen's slope but no intercept...?
slope_method = "SS"
if slope_method == "SS": slopes = SS_slopes
elif slope_method == "LR": slopes = LR_slopes
else: print("error with slope method")
summary_stats_all = do_summary_stats(slopes, slope_method, dat_ts, dat_davg, dat_numerictime, dat_seas, dat_timenumeric_seas)

# ====  Parametric Trend Detection =====
# use LR or SS with adjusted ESS from DST
type_CI = "parametric"
DTS_CI_out = []
alpha = 0.05
if type_CI == "parametric":
    DTS_CI_out = dst_trend_sig(summary_stats_all, dat_ts, alpha, dat_davg, dat_seas, dat_numerictime, dat_timenumeric_seas)

# ==== Non-Parametric Trend Detection ====
# - MK w/ Prewhitening (e.g., 3PW)
type_CI = "nonparametric"
resolution = 0.0001  # of measurements to use in ties calc in Kendalls Tau and S
if type_CI == "nonparametric":
    print("MK w/ Sens and 3PW method")
    mk_3pw_out = mk3pw_trend_sig(seasons_per_year, resolution, dat_davg, dat_pytime, dat_seas, dat_pytime_seas)

print(mk_3pw_out)





# for each season individually, and all as a whole, do:
# Detrend
# Test residuals for normality
# Test residuals for heteroskedacity
# use FFT to test resids for periocity
# compute ACF and AR1 Coeff

# ==== Test Normality ====
# Test resid for norm using Shapiro-Wilks
#residuals = dat_davg_dtrnd
# stat, p_value = stats.shapiro(residuals)
# print(f'Shaprio Wilks Statistic: {stat}, p-value: {p_value}')
# if p_value > 0.05:
#     print("Residuals likely are normally distributed. Parametric methods could be used for CI estimation.")
#     type_CI = "parametric"
# else:
#     print("Residuals are unlikely to be drawn from a normal distribution. Non-parametric methods for CI estimation should be considered.")
#     type_CI = "non-parametric"

# ==== Tests for Heteroscedacity ====
# White's test
# _, p_value_whites, _, _ = het_white(residuals, exog=add_constant(dat_numerictime))
# print(f"White's test p-value (low = unlikely heteroscedastic): {p_value_whites}")
# # Goldfeld-Quandt test
# _, p_value_goldfield, _ = het_goldfeldquandt(residuals, add_constant(dat_numerictime))
# print(f'Goldfeld-Quandt test p-value (low = unlikely heteroscedastic): {p_value_goldfield}')

# ==== Test for Seasonality using FFT ====
# freq, power_spec = do_fft(dat_davg_dtrnd, len(dat_davg_dtrnd))
# peak_indices = np.argsort(power_spec)[::-1][:5]  # Top 5 peaks
# peak_frequencies = np.abs(freq[peak_indices])    # corresponding frequencies
# peak_freq_period_yrs = 1/peak_frequencies*1/dat_ts

# plot power spec
# len_ts = len(dat_davg_dtrnd)
# plt.plot(freq[:len_ts // 2], power_spec[:len_ts // 2]) # Exclude the negative freqs
# plt.title('Power Spectrum')
# plt.xlabel('Frequency (cycles per data point)')
# plt.ylabel('Power')
# plt.show()

# check for sub-annual seasonal periocity
# print("Top peak frequencies in years:", peak_freq_period_yrs)
# print(np.any(peak_freq_period_yrs <= 1))
# subannual_freq_tf = np.any(peak_freq_period_yrs <= 1)
# if subannual_freq_tf:
#     print("Subannual cycles are detected in anomaly residuals.")
# else:
#     print("Subannual cycles are not detected in anomaly residuals.")

# ==== Compute Autocorrelation Function =====
# alpha = 0.05
# nlags = min(len(residuals)-1, 96)
# acf_values = acf(residuals.values, adjusted=True, fft=False, alpha=alpha, nlags=nlags, missing="conservative")
# plot_acf(residuals.data, adjusted=True, fft=False, alpha=alpha, lags=nlags)
# plt.title('Autocorrelation Function with Confidence Intervals')
# plt.xlabel('Lag')
# plt.ylabel('Autocorrelation (-1, 1)')
# plt.show()
# AR1_coef = acf_values[0][1]
# print("The AR(1) Coefficient:", AR1_coef)

    # #
    # # print("proceeding with CI estimation using decorrelation time scale")
    # DTS = (1 + AR1_coef) / (1 - AR1_coef) # decorr. time scale in ts units
    # ESS = len(dat_davg) / DTS             # effective sample size
    #
    # print("The Decorrelation Time Scale (DTS; units: yr): ", DTS / dat_ts)
    # print("The Effective Sample Size (ESS): ", ESS)
    #
    # for season in range(0, seasons_per_year):
    #     dat_davg_np = np.asarray(dat_seas[season].values)
    #     d1 = seasonal_dates_numeric[season]
    #     slope_dst, inter_dst, r_val_dst, p_val_dst, std_err_dst, _ = linregress_GO(d1,
    #                                                                                dat_davg_np, ESS)
    #     df = ESS - 2
    #
    #     t_critical = stats.t.ppf((1 + (1 - alpha)) / 2, df)
    #     margin_of_error = t_critical * std_err_dst
    #     confidence_interval = ((slope_dst - margin_of_error) * 365, (slope_dst + margin_of_error) * 365)
    #     slope_dst = slope_dst * 365
    #     DTS_CI_out.append({"season " + str(season):{"slope": slope_dst, "inter": inter_dst,
    #                                                 "r_val": r_val_dst, "p_val": p_val_dst,
    #                                                 "std_err": std_err_dst, "slope_lcl": confidence_interval[0],
    #                                                 "slope_ucl": confidence_interval[1]}})
    #     print("Slope from LR for season ", str(season), ": ", str(slope_dst), " lcl:", str(confidence_interval[0]),
    #           " ucl:", str(confidence_interval[1]), " pval (DST):", str(p_val_dst))
    # # adjust sample size downward
    # slope_dst, inter_dst, r_val_dst, p_val_dst, std_err_dst, _ = linregress_GO(dat_numerictime.values,
    #                                                                            dat_davg.values,
    #                                                                            ESS)
    # df = ESS - 2
    # t_critical = stats.t.ppf((1 + (1-alpha)) / 2, df)
    # margin_of_error = t_critical * std_err_dst
    # confidence_interval = ((slope_dst - margin_of_error) * 365, (slope_dst + margin_of_error) * 365)
    # slope_dst = slope_dst * 365
    # DTS_CI_out.append({"all seasons": {"slope": slope_dst, "inter": inter_dst,
    #                                    "r_val": r_val_dst, "p_val": p_val_dst,
    #                                    "std_err": std_err_dst, "slope_lcl": confidence_interval[0],
    #                                    "slope_ucl": confidence_interval[1]}})
    # print(f"Slope: {slope_dst*365}")
    # print(f"Confidence Interval (95%): {confidence_interval}")
    # print(f"P val:  {p_val_dst}")



    # if seasons_per_year == 1:
    #     out = mk.mk_temp_aggr(dat_pytime, np.asarray(dat_davg), resolution)
    #     mk_3pw_out.append({"Season " + str(season): out})
    #     print(out)
    # else:
    #     print("Seasons analysed individually:")
    #     for season in range(1, seasons_per_year + 1):
    #         out = mk.mk_temp_aggr([seasonal_dates[season-1]], [dat_seas[season-1]], resolution)
    #         mk_3pw_out.append({"season " + str(season): out[0]})
    #         print(season, ": ", out)
    #     out = mk.mk_temp_aggr(seasonal_dates, dat_seas, resolution)
    #     mk_3pw_out.append({"all seasons": out})
    #     print("Seasons taken together:")
    #     for n in range(seasons_per_year):
    #         print('Season {ind}:'.format(ind=n + 1), out[n])
    #     print('Overall: ', out[len(out)-1])

#return: all seasons together, separately - slope, upper, lower cl, p val, ar1 coeff
# to do, bootstrap
# else - to do