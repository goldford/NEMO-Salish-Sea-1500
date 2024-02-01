# misc tools by GO for the trend analysis
# by GO Jan 2024

import numpy as np
import netCDF4 as nc
import os
from scipy.fft import fft, ifft
import xarray as xr
import pandas as pd
import mannkendall as mk
from scipy import stats
from scipy.stats import linregress
from linregress_GO import linregress_GO
from statsmodels.stats.diagnostic import het_goldfeldquandt, het_white
from statsmodels.tools.tools import add_constant
from statsmodels.tsa.stattools import acf

def mk3pw_trend_sig(seasons_per_year, resolution, alpha_mk, dat, dat_pytime, dat_seas, dat_pytime_seas):
    #created Jan 2024 by G Oldford
    # purpose: implement the '3PW' approach to prewhitening to account for autocorrelation
    #          as described in Collaud Coen et al 2020; https://doi.org/10.5194/amt-13-6945-2020
    # inputs:
    #    seasons_per_year  - integer; number of seasons in a year
    #    resolution        - float; the precision of the measurements in the data
    #    alpha_mk          - the confidence level threshold for mann-kendall (e.g., 0.95 for 95%)
    #    dat               - DataArray; xarray dataarray with one var
    #    dat_pytime        - ndarray; numpy array of dates as python datetime, as required by mk
    #    dat_seas          - DataArray; list of xarray DataArrays with seasonal data
    #    dat_pytime_seas   - list; list of ndarray numpy arrays with python datetimes, one list element for each season
    # returns:
    #    mk_3pw_out  - a list of dictionary elements with stats outputs

    mk_3pw_out = []
    print("MK w/ Sens and 3PW method")
    if seasons_per_year == 1:
        out = mk.mk_temp_aggr(dat_pytime, np.asarray(dat), resolution, alpha_mk=alpha_mk)
        mk_3pw_out.append({"season " + str(1): out})
    else:
        print("Seasons analysed individually:")
        for season in range(1, seasons_per_year + 1):
            out = mk.mk_temp_aggr([dat_pytime_seas[season - 1]], [dat_seas[season - 1]],
                                  resolution, alpha_mk=alpha_mk)
            mk_3pw_out.append({"season " + str(season): out[0]})
        out = mk.mk_temp_aggr(dat_pytime_seas, dat_seas, resolution)
        mk_3pw_out.append({"all seasons": out})
    return mk_3pw_out

def dst_trend_sig(summary_stats_all, dat_ts, alpha, dat, dat_seas, dat_numerictime, dat_numerictime_seas):
    # created by G Oldford Jan 29
    # purpose: use decorrelation time scale to estimate confidence intervals around secular trend slope
    # inputs:
    #    summary_stats_all    - list; list of dictionary elements with stats necessary to run the trend analysis
    #    dat_ts               - integer; data time scale, the number of data points within a year (i.e., 'seasons')
    #    alpha                - float; threshold for conf. level using t-test (e.g., 0.05 for 95% CL, two tailed)
    #    dat                  - DataArray; xarray dataarray with one var
    #    dat_seas             - list; list of xarray data arrays, one list element for each season
    #    dat_numerictime      - DataArray; xarray datarray of time as numeric units (days)
    #    dat_numerictime_seas - list; list of xarray dataarrays, one list element for each season
    # returns:
    #    DTS_CI_out  - a list of dictionary elements with results of the decorrelation time scale analysis of trend sig

    DTS_CI_out = []
    for seas in range(0, len(summary_stats_all)):
        for key in summary_stats_all[seas].keys():
            AR1_coef = summary_stats_all[seas][key]['acf-ar1-coeff']
            DTS = (1 + AR1_coef) / (1 - AR1_coef)  # decorr. time scale in ts units
            if key == 'all seasons':
                ESS = len(dat) / DTS  # effective sample size
                dat_np = np.asarray(dat.values)
                dat_ntime = dat_numerictime
                DTS_yr = DTS / dat_ts

            else:
                i = 0
                for s2 in dat_seas:
                    if i == seas:
                        ESS = len(dat_seas[i]) / DTS
                        dat_np = np.asarray(dat_seas[i].values)
                        dat_ntime = dat_numerictime_seas[i]
                        DTS_yr = DTS
                    i += 1
            print("The Decorrelation Time Scale (DTS; units: yr): ", DTS_yr)
            print("The Effective Sample Size (ESS): ", ESS)

            slope_dst, inter_dst, r_val_dst, p_val_dst, std_err_dst, _ = linregress_GO(dat_ntime, dat_np, ESS)
            df = ESS - 2
            t_critical = stats.t.ppf((1 + (1 - alpha)) / 2, df)
            margin_of_error = t_critical * std_err_dst
            confidence_interval = ((slope_dst - margin_of_error) * 365, (slope_dst + margin_of_error) * 365)
            slope_dst = slope_dst * 365
            DTS_CI_out.append({key: {"slope": slope_dst, "inter": inter_dst,
                                     "r_val": r_val_dst, "p_val": p_val_dst,
                                     "std_err": std_err_dst, "DTS_yr": DTS_yr,
                                     "ESS": ESS,
                                     "slope_lcl": confidence_interval[0],
                                     "slope_ucl": confidence_interval[1]}})
            print("Slope from LR for season ", key, ": ", str(slope_dst), " lcl:", str(confidence_interval[0]),
                  " ucl:", str(confidence_interval[1]), " pval (DST):", str(p_val_dst))

    return DTS_CI_out

def do_summary_stats(slopes, slope_method, dat_ts, dat, dat_timenumeric, dat_seas, dat_timenumeric_seas,
                     alpha_acf=0.05):
    # created by G Oldford Jan 2024
    # purpose: generate a dictionary of summary stats using detrended anomalies (seasonal and annual)
    #          meant to help with understanding how best to test for trend and significance
    # input:
    #     slopes               - list; list of dictionary elements containing slopes
    #     slope_method         - string; "SS" (Theil-Sen) or "LR" (linear regression)
    #     dat_ts               - integer; data time scale - number of seasons per year (24, 12, 4, 1)
    #     dat                  - dataArray; xr DataArray of data, annual, with one variable
    #     dat_timenumeric      - dataArray; list of dates, annual, timedelta
    #     dat_seas             - list; list of xr DataArrays for each season
    #     dat_timenumeric_seas - list; list of xr Datarrays with numeric dates for each season, timedelta
    # returns:
    #     summary_stats_all    - dictionary with a slew of exploratory test results

    slope = -999
    res = -999
    summary_stats_all = []
    for seas in range(0, len(slopes)):
        print(slopes[seas])
        print(slopes[seas].keys())

        keys = slopes[seas].keys()
        for key in slopes[seas].keys():
            if key == 'all seasons':
                slope = slopes[seas]['all seasons']['slope']
                inter = slopes[seas]['all seasons']['inter']
                dat_1seas = dat
                dat_1seas_timenum = dat_timenumeric
            else:
                slope = slopes[seas]['season ' + str(seas + 1)]['slope']
                inter = slopes[seas]['season ' + str(seas + 1)]['inter']
                dat_1seas = dat_seas[seas]
                dat_1seas_timenum = dat_timenumeric_seas[seas]

            if slope_method == "SS":
                residuals = dat_1seas - (slope * (dat_timenumeric / 365) + inter)
            # else
            # to do - different time and time mult if lr

            stat_sw, p_value_sw = stats.shapiro(residuals)  # test for normality (small p val = small prob of normal)
            _, p_value_whites, _, _ = het_white(residuals, exog=add_constant(dat_1seas_timenum))  # test for skew
            _, p_value_goldfield, _ = het_goldfeldquandt(residuals, add_constant(dat_1seas_timenum))  # test for skew
            # test periodicity using fft
            freq, power_spec = do_fft(residuals, len(residuals))
            peak_indices = np.argsort(power_spec)[::-1][:5]  # Top 5 peaks
            peak_frequencies = np.abs(freq[peak_indices])  # corresponding frequencies

            if key == "all seasons":
                peak_freq_period_yrs = 1 / peak_frequencies * 1 / dat_ts
            else:
                peak_freq_period_yrs = 1 / peak_frequencies

            # autocorrelation
            nlags = len(residuals - 1)
            acf_values = acf(residuals.values, adjusted=True, fft=False, alpha=alpha_acf, nlags=nlags,
                             missing="conservative")
            AR1_coef = acf_values[0][1]

            summary_stats_all.append({key: {"norm-shapiro-wilks stat": stat_sw,
                                    "norm-shapiro-wilks pval": p_value_sw,
                                    "het-white pval": p_value_whites,
                                    "het-goldfield pval": p_value_goldfield,
                                    "fft-top5-peak-freqs (yrs)": peak_freq_period_yrs,
                                    "fft-freq": freq,
                                    "fft-power-spec": power_spec,
                                    "acf-ar1-coeff": AR1_coef,
                                    "acf-all-nlag-coeffs": acf_values
                                    }
                              })
    return summary_stats_all

def prep_data_seas(dat, dat_pytime, dat_numerictime, seasons_per_yr, time_nm):
    # created by G Oldford Jan 2024
    # does lin regression for data grouped by season and annually
    # inputs
    #   dat             - depth averaged data in xarray with one variable
    #   dat_pytime      - time as python datetime
    #   dat_numerictime - time, numeric, as timedelta e.g. dat_davg[time_nm] - start_date) / np.timedelta64(1, 'D')
    #   seasons_per_yr  - integer, number of seasons per year
    #   time_nm         - the dimension name corresponding to time
    # returns
    #   LR_slopes       - a list of dictionary elements

    seasonal_dates = []
    seasonal_dates_numeric = []
    dat_seas = []  # time series as list
    if seasons_per_yr == 1:
        seasonal_dates.append((dat_pytime))
        seasonal_dates_numeric.append((dat_numerictime))
        dat_seas.append((dat))
    else:
        for season in range(1, seasons_per_yr + 1):
            if seasons_per_yr == 24:
                indices = np.where(
                    (dat[time_nm].dt.dayofyear >= (season * 15) - 10) &
                    (dat[time_nm].dt.dayofyear < ((season + 1) * 15) - 10))[0]
            elif seasons_per_yr == 12:
                indices = np.where(dat[time_nm].dt.month == season)
            elif seasons_per_yr == 4:
                if season == 1:
                    season_name = "winter"
                elif season == 2:
                    season_name = "spring"
                elif season == 3:
                    season_name = "summer"
                elif season == 4:
                    season_name = "fall"
                else:
                    print("error w/ seasons")
                indices = np.where(dat['season'] == season_name)
            seasonal_dates.append((dat_pytime[indices]))
            seasonal_dates_numeric.append((dat_numerictime[indices]))
            dat_seas.append((dat[indices]))

    return dat_seas, seasonal_dates, seasonal_dates_numeric

def do_lr_seasons(dat, dat_seas, dat_time, dat_time_seas, seas_per_yr):
    # created by G Oldford Jan 2024
    # does lin regression for data grouped by season and annually
    # inputs
    #   dat           - depth averaged data in xarray with one variable
    #   dat_seas      - same as above but a list of the data by season
    #   dat_time      - time, numeric, as timedelta e.g. dat_davg[time_nm] - start_date) / np.timedelta64(1, 'D')
    #   dat_time_seas - same as above but seasonal
    # returns
    #   LR_slopes     - a list of dictionary elements
    LR_slopes = []
    for season in range(0, seas_per_yr):
        d1 = dat_time_seas[season]
        dmult = dat_seas[season]
        slope, inter, rval, pval, stderr = linregress(d1, dmult.values)
        LR_slopes.append({"season "+str(season+1): {"slope": slope * 365, "inter": inter * 365,
                          "rval": rval, "pval": pval, "stderr": stderr}})

    slope, inter, rval, pval, _ = linregress(dat_time.values, dat.values)
    slope = slope * 365  # returns slope in /day (the units of time_numeric
    inter = inter * 365  # returns slope in /day (the units of time_numeric
    LR_slopes.append({"all seasons": {"slope": slope, "inter": inter,
                                      "rval": rval, "pval": pval, "stderr": stderr}})
    return LR_slopes

def do_ss_seasons(dat, dat_seas, dat_pytime, dat_pytime_seas, seasons_per_yr, resolution):
    # created by G Oldford Jan 2024
    # calling mks.sens_slope is not usually called directly, so this code was borrowed from mannkendall lib
    # Intercept calculated using Conover, W.J. (1980) method (from pymannkendall)
    #
    # inputs
    #   dat              - depth-averaged xarray with one variable and time dimension
    #   dat_seas         - a version of the above, but a list of d-avg values by season
    #   dat_pytime       - xarray or numpy array, the time counter var from dat_davg as py datetime, required by mk
    #   dat_pytime_seas  - a version of the above, but a list of time vals by season
    #   seasons_per_yr   - integer, number of seasons per year
    #   resolution       - the precision of the measurements, significant digits
    # returns
    #   list of dictionary elements, slope and upper conf. limit, lower conf, intercept - both annual and seasonal
    #
    #  note all slopes are converted to units of years but sen_slope returns units of seconds

    SS_slopes = []
    alpha_cl = 90

    for season in range(0, seasons_per_yr):
        dat_davg_np = np.asarray(dat_seas[season])
        t = mk.mkt.nb_tie(dat_davg_np, resolution)
        (s, n) = mk.mks.s_test(dat_davg_np, dat_pytime_seas[season])
        k_var = mk.mkt.kendall_var(dat_davg_np, t, n)
        z = mk.mks.std_normal_var(s, k_var)
        mk_out = mk.mks.sen_slope(dat_pytime_seas[season], dat_davg_np, k_var, alpha_cl=alpha_cl)  # units of second!
        mult = 3600 * 24 * 365.25
        slope_yr = mk_out[0] * mult
        slope_yr_min = mk_out[1] * mult
        slope_yr_max = mk_out[2] * mult
        inter_yr = np.nanmedian(dat_davg_np) - np.median(
            np.arange(len(dat_davg_np))[~np.isnan(dat_davg_np.flatten())]) * slope_yr
        SS_slopes.append({"season " + str(season+1): {"slope": slope_yr,
                                                      "slope_LCL": slope_yr_min,
                                                      "slope_UCL": slope_yr_max,
                                                      "inter":inter_yr}})
    dat_davg_np = np.asarray(dat)
    t = mk.mkt.nb_tie(dat_davg_np, resolution)
    (s, n) = mk.mks.s_test(dat_davg_np, dat_pytime)
    k_var = mk.mkt.kendall_var(dat_davg_np, t, n)
    z = mk.mks.std_normal_var(s, k_var)
    mk_out = mk.mks.sen_slope(dat_pytime, dat_davg_np,k_var,alpha_cl=alpha_cl) # units of second!
    mult = 3600 * 24 * 365.25
    slope_yr = mk_out[0] * mult; slope_yr_min = mk_out[1] * mult; slope_yr_max = mk_out[2] * mult
    inter_yr = np.nanmedian(dat_davg_np) - np.median(np.arange(len(dat_davg_np))[~np.isnan(dat_davg_np.flatten())]) * slope_yr
    # or median(x) - (n-1)/2 *slope
    SS_slopes.append({"all seasons": {"slope": slope_yr,
                                      "slope_LCL": slope_yr_min,
                                      "slope_UCL": slope_yr_max,
                                      "inter":inter_yr}})
    return SS_slopes

def do_seasonal_avg(dat, time_inc, time_dim="time_counter"):
    # created by G Oldford Jan 2024
    # purpose: creates second seasonal version of the original array
    #

    if time_inc == "annual":
        time_dim_grp = time_dim + ".year"
        dat = dat.groupby(time_dim_grp).mean(dim=time_dim)
        date_range = [np.datetime64(f'{year}-01-01') for year in dat['year'].values]
        dat['year'] = date_range
        dat = dat.rename({'year': time_dim})
        time_nm = time_dim
        dat_ts = 1
    elif time_inc == "monthly":
        dat[time_dim] = dat[time_dim].dt.strftime('%Y-%m')
        dat = dat.groupby(time_dim).mean(dim=time_dim)
        date_range = pd.to_datetime(dat[time_dim], format='%Y-%m')
        dat[time_dim] = date_range
        dat_ts = 12
        time_nm = time_dim
    elif time_inc == "seasonal":
        # awkward but no better way found
        time_dim_grp_mo = time_dim + ".month"
        time_dim_grp_yr = time_dim + ".year"
        seasons = {'winter': [12, 1, 2], 'spring': [3, 4, 5], 'summer': [6, 7, 8], 'fall': [9, 10, 11]}
        dat['time_counter_seas'] = dat[time_dim].values[0]
        start_yr = dat[time_dim].dt.year.values[0]
        end_yr = dat[time_dim].dt.year.values[-1]
        for yr in range(start_yr, end_yr + 1, 1):
            for season_name, season_months in seasons.items():
                date_value = np.datetime64(f'{yr}-{season_months[1]:02}-01')

                dat['time_counter_seas'] = xr.where(
                    (dat[time_dim_grp_mo].isin(season_months) & (dat[time_dim_grp_yr] == yr)),
                    date_value,
                    dat['time_counter_seas'])

        dat = dat.groupby('time_counter_seas').mean(dim=time_dim)

        # second loop to add a 'season' label for later
        dat['season'] = 'other'  # initialise
        for yr in range(start_yr, end_yr, 1):
            for season_name, season_months in seasons.items():
                dat['season'] = xr.where(
                    dat['time_counter_seas.month'].isin(season_months),
                    season_name,
                    dat['season'])

        dat_ts = 4
        time_nm = 'time_counter_seas'
    elif time_inc == 'biweekly':
        dat_ts = 24
        time_nm = time_dim
    else:
        print('Unrecognized time increment')

    return dat, dat_ts, time_nm

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