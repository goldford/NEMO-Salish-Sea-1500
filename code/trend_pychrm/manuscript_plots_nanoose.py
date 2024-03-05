# Created Feb 2024 by G Oldford
# Purpose: produce plots of model vs obs at Nanoose stn

from trends_test_GO import do_prelim_trend_analysis
from GO_tools import get_dat, write_mk3pw_data_to_csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import numpy as np

# ==== Paths ====
meshm_p = '../../data/mesh mask/'
meshm_f = 'mesh_mask_20210406.nc'
dat_p = '../climatol_intermediate_files/'

# this does nothing yet except select files, imperfectly
year_min = 1980
year_max = 2018

#dat_f = 'Nanoose_obs_anom_temp_' + str(year_min) + '-' + str(year_max) + '.nc'
#dat_f = 'Nanoose_obs_anom_temp_' + str(year_min) + '-' + str(year_max) + '.nc'
dat_f_obs = 'Nanoose_obs_anom_temp_' + str(year_min) + '-' + str(year_max) + '.nc'
dat_f_obs_full = 'Nanoose_obs_anom_temp_' + str(1970) + '-' + str(year_max) + '.nc'
dat_f_mod1 = 'RUN203mod_anom_temp_1980-2018.nc'
dat_f_mod2 = 'RUN216mod_anom_temp_1980-2018.nc'

# ==== Params for Analysis ====
depth_min = 4.5
depth_max = 400
var = 'temperature'
remove_nans = True
use_abs = False # make sure false if you want trends in variable (not variability)!!
alpha_DTS_CI = 0.05
alpha_MK = 95
time_inc = 'seasonal' # biweekly, monthly, seasonal, annual (manuscipr seasonal)
resolution = 0.0001  # of measurements to use in ties calc in Kendalls Tau and S
if time_inc == 'annual': seasons_per_year = 1
elif time_inc == 'seasonal': seasons_per_year = 4
elif time_inc == 'monthly': seasons_per_year = 12
elif time_inc == 'biweekly': seasons_per_year = 24

# get dat does depth averaging
(d_o, d_pt_o, d_nt_o,
 d_se_o, d_pt_se_o, d_tn_se_o,
 d_ts, time_dim, se_year) = get_dat(meshm_p, meshm_f,
                                    dat_p, dat_f_obs, var, time_inc, use_abs,
                                    dep_int=True, depth_min=depth_min, depth_max=depth_max)

(d_o_full, d_pt_o_full, d_nt_o_full,
 d_se_o_full, d_pt_se_o_full, d_tn_se_o_full,
 d_ts, time_dim, se_year) = get_dat(meshm_p, meshm_f, dat_p, dat_f_obs_full, var, time_inc, use_abs,
                                    dep_int=True, depth_min=depth_min, depth_max=depth_max)

(d_m1, d_pt_m1, d_nt_m1,
 d_se_m1, d_pt_se_m1, d_tn_se_m1,
 d_ts, time_dim, se_year) = get_dat(meshm_p, meshm_f, dat_p, dat_f_mod1, var, time_inc, use_abs,
                                    dep_int=True, depth_min=depth_min, depth_max=depth_max)

(d_m2, d_pt_m2, d_nt_m2,
d_se_m2, d_pt_se_m2, d_tn_se_m2,
 d_ts, time_dim, se_year) = get_dat(meshm_p, meshm_f, dat_p, dat_f_mod2, var, time_inc, use_abs,
                                    dep_int=True, depth_min=depth_min, depth_max=depth_max)

# ==== Deal with missing values ====
if remove_nans:
    nan_mask = np.logical_or(np.logical_or(np.isnan(d_o), np.isnan(d_m1)), np.isnan(d_m2))
    d_o[nan_mask] = np.nan; d_m1[nan_mask] = np.nan; d_m2[nan_mask] = np.nan
    #d_pt_o[nan_mask] = np.nan; d_pt_m1[nan_mask] = np.nan; d_pt_m2[nan_mask] = np.nan
    #d_nt_o[nan_mask] = np.nan; d_nt_m1[nan_mask] = np.nan; d_nt_m2[nan_mask] = np.nan
    for s in range(0,len(d_se_o)-1):
        nan_mask = np.logical_or(np.logical_or(np.isnan(d_se_o[s]), np.isnan(d_se_m1[s])), np.isnan(d_se_m2[s]))
        d_se_o[s][nan_mask] = np.nan; d_se_m1[s][nan_mask] = np.nan; d_se_m2[s][nan_mask] = np.nan
        #d_pt_se_o[s][nan_mask] = np.nan; d_pt_se_m1[s][nan_mask] = np.nan; d_pt_se_m2[s][nan_mask] = np.nan
        #d_tn_se_o[s][nan_mask] = np.nan; d_tn_se_m1[s][nan_mask] = np.nan; d_tn_se_m2[s][nan_mask] = np.nan

    # this returns depth averaged data and time arrays, with and without seasonal arrangement
# and returns the slopes, summary stats, and sig tests using two methods

(table_obs, summary_stats_all_obs,
 mk_3pw_out_obs, DTS_CI_out_obs,
 LR_slopes_obs, SS_slopes_obs) = do_prelim_trend_analysis(d_o, d_pt_o, d_nt_o,
                                                          d_se_o, d_pt_se_o, d_tn_se_o,
                                                          se_year, d_ts, resolution,
                                                          alpha_DTS_CI, alpha_MK)

(table_obs_full, summary_stats_all_obs_full,
 mk_3pw_out_obs_full, DTS_CI_out_obs_full,
 LR_slopes_obs_full, SS_slopes_obs_full) = do_prelim_trend_analysis(d_o_full, d_pt_o_full,
                                                                    d_nt_o_full, d_se_o_full,
                                                                    d_pt_se_o_full, d_tn_se_o_full,
                                                                    se_year, d_ts, resolution,
                                                                    alpha_DTS_CI, alpha_MK)

(table_mod1, summary_stats_all_mod1,
 mk_3pw_out_mod1, DTS_CI_out_mod1,
 LR_slopes_mod1, SS_slopes_mod1) = do_prelim_trend_analysis(d_m1, d_pt_m1, d_nt_m1,
                                                            d_se_m1, d_pt_se_m1, d_tn_se_m1,
                                                            se_year, d_ts, resolution,
                                                            alpha_DTS_CI, alpha_MK)

(table_mod2, summary_stats_all_mod2,
 mk_3pw_out_mod2, DTS_CI_out_mod2,
 LR_slopes_mod2, SS_slopes_mod2) = do_prelim_trend_analysis(d_m2, d_pt_m2, d_nt_m2,
                                                            d_se_m2, d_pt_se_m2, d_tn_se_m2,
                                                            se_year, d_ts, resolution,
                                                            alpha_DTS_CI, alpha_MK)

if use_abs:abs_lab = 'absolute'
else: abs_lab = ''
write_mk3pw_data_to_csv(mk_3pw_out_obs, 'NanooseObs_1980-2018', 'mk_3pw_out_' + str(depth_min) + 'to' + str(depth_max) + time_inc + '_Nanoose1980-2018' + abs_lab + '.csv')
write_mk3pw_data_to_csv(mk_3pw_out_obs_full, 'NanooseObs_1970-2018', 'mk_3pw_out_' + str(depth_min) + 'to' + str(depth_max) + time_inc + '_Nanoose1970-2018' + abs_lab + '.csv')
write_mk3pw_data_to_csv(mk_3pw_out_mod1, 'HOTSS_v1.01', 'mk_3pw_out_' + str(depth_min) + 'to' + str(depth_max) + time_inc + '_HOTSSv101' + abs_lab + '.csv')
write_mk3pw_data_to_csv(mk_3pw_out_mod2, 'HOTSS_v1.02', 'mk_3pw_out_' + str(depth_min) + 'to' + str(depth_max) + time_inc + '_HOTSSv102' + abs_lab + '.csv')

# plot trends by season
if time_inc == 'seasonal':
    import pandas as pd
    fig, ax = plt.subplots(figsize=(5,2))
    stats = SS_slopes_obs
    means = [stats[0]['season 1']['slope']*10,
                   stats[1]['season 2']['slope']*10,
                   stats[2]['season 3']['slope']*10,
                   stats[3]['season 4']['slope']*10]
    upper_limit = [stats[0]['season 1']['ucl']*10,
                   stats[1]['season 2']['ucl']*10,
                   stats[2]['season 3']['ucl']*10,
                   stats[3]['season 4']['ucl']*10]
    lower_limit = [stats[0]['season 1']['lcl']*10,
                   stats[1]['season 2']['lcl']*10,
                   stats[2]['season 3']['lcl']*10,
                   stats[3]['season 4']['lcl']*10]
    seasons = np.asarray(['Winter','Spring','Summer','Fall'])

    lcl = np.asarray(lower_limit)
    ucl = np.asarray(upper_limit)
    means = np.asarray(means)
    err_lcl = np.abs(means - lcl)
    err_ucl = np.abs(ucl - means)

    # Plotting
    seasons_idx = range(len(seasons))
    seasons_idx = [i - 0.1 for i in seasons_idx]  # offset
    ax.plot(means, seasons_idx, marker='o', linestyle='', color='r', markerfacecolor='w')
    ax.plot([lcl[0],ucl[0]],[seasons_idx[0],seasons_idx[0]], color='r', label='Obs.')
    ax.plot([lcl[1], ucl[1]], [seasons_idx[1], seasons_idx[1]], color='r')
    ax.plot([lcl[2], ucl[2]], [seasons_idx[2], seasons_idx[2]], color='r')
    ax.plot([lcl[3], ucl[3]], [seasons_idx[3], seasons_idx[3]], color='r')

    stats = SS_slopes_mod2
    means = [stats[0]['season 1']['slope'] * 10,
             stats[1]['season 2']['slope'] * 10,
             stats[2]['season 3']['slope'] * 10,
             stats[3]['season 4']['slope'] * 10]
    upper_limit = [stats[0]['season 1']['ucl'] * 10,
                   stats[1]['season 2']['ucl'] * 10,
                   stats[2]['season 3']['ucl'] * 10,
                   stats[3]['season 4']['ucl'] * 10]
    lower_limit = [stats[0]['season 1']['lcl'] * 10,
                   stats[1]['season 2']['lcl'] * 10,
                   stats[2]['season 3']['lcl'] * 10,
                   stats[3]['season 4']['lcl'] * 10]
    seasons = np.asarray(['Winter', 'Spring', 'Summer', 'Fall'])

    lcl = np.asarray(lower_limit)
    ucl = np.asarray(upper_limit)
    means = np.asarray(means)
    err_lcl = np.abs(means - lcl)
    err_ucl = np.abs(ucl - means)

    # Plotting
    seasons_idx = range(len(seasons))
    seasons_idx = [i + 0.1 for i in seasons_idx] # offset
    ax.plot(means, seasons_idx, marker='o', linestyle='', color='g', markerfacecolor='w')
    ax.plot([lcl[0], ucl[0]], [seasons_idx[0], seasons_idx[0]], linestyle='--', color='g', label='Mod.')
    ax.plot([lcl[1], ucl[1]], [seasons_idx[1], seasons_idx[1]], linestyle='--', color='g')
    ax.plot([lcl[2], ucl[2]], [seasons_idx[2], seasons_idx[2]], linestyle='--', color='g')
    ax.plot([lcl[3], ucl[3]], [seasons_idx[3], seasons_idx[3]], linestyle='--', color='g')

    ax.set_xlim(-0.1,0.2)
    ax.set_yticks(seasons_idx)
    ax.set_yticklabels(['Winter','Spring','Summer','Fall'])
    ax.set_xlabel('Temp. Anom. Trend (\u00B0C) / decade)')

    # Customize y-axis labels
    # err_lcl = means - list(lower_limit.values())
    # err_ucl = list(upper_limit.values()) - means
    # err_lcl = np.abs(err_lcl)
    # err_ucl = np.abs(err_ucl)
    # ax.errorbar(seasons, means, yerr=[err_lcl,err_ucl], fmt='')
    #ax.plot(seasons, means, marker='o', linestyle='', color='blue')
    # Customize x-axis labels
    # plt.xticks(range(len(seasons)), seasons)
    ax.invert_yaxis()
    if use_abs == True: abs_lab='_abs'
    else: abs_lab = ''
    plt.savefig('modobsnanoose_' + str(depth_min) + '-' + str(depth_max) + '_' + abs_lab + '_' + time_inc + '.png', dpi=300)
    plt.legend()
    plt.show()

fig, ax = plt.subplots()
ax.fill_between(d_pt_o, d_o, d_m1, color='grey', label='', alpha=0.7, zorder=0)
#ax.scatter(dat_pytime, dat_davg, label='Obs', s=1, color='r', marker='o')
#ax.scatter(dat_pytime, mod_davg, label='Mod', s=1, color='g', marker='x')
ax.plot(d_pt_o, d_o, label='Obs', color='k', linestyle='-', linewidth=1)
ax.plot(d_pt_m1, d_m1, label='Mod', color='r', linestyle='--', linewidth=1)
ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.set_ylim(-1,1)

plt.title("HOTSS v1.01 Modelled vs. Observed Temperature Anomalies at Nanoose Stn, " + str(depth_min) + " to " + str(depth_max) + " m")
plt.xlabel('Time')
plt.ylabel('Anomaly (deg C)')
plt.legend()
plt.show()

fig, ax = plt.subplots()
ax.fill_between(d_pt_o, d_o, d_m2, color='grey', label='', alpha=0.7, zorder=0)
#ax.scatter(dat_pytime, dat_davg, label='Obs', s=1, color='r', marker='o')
#ax.scatter(dat_pytime, mod_davg, label='Mod', s=1, color='g', marker='x')
ax.plot(d_pt_o, d_o, label='Obs', color='k', linestyle='-', linewidth=0.8)
ax.plot(d_pt_m2, d_m2, label='Mod', color='r', linestyle='--', linewidth=0.8)
ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.set_ylim(-1,1)
plt.title("HOTSS v1.02 Modelled vs. Observed Temperature Anomalies at Nanoose Stn, " + str(depth_min) + " to " + str(depth_max) + " m")
plt.xlabel('Time')
plt.ylabel('Anomaly (deg C)')
plt.legend()
plt.show()

fig, ax = plt.subplots(figsize=(7.5,4.5))

# optional bit - using moving average for tiny publication plots
# window_size = 2
# # Apply the moving average using numpy.convolve
# d_o = np.convolve(d_o, np.ones(window_size)/window_size, mode='valid')
# d_o_full = np.convolve(d_o_full, np.ones(window_size)/window_size, mode='valid')
# d_m2 = np.convolve(d_m2, np.ones(window_size)/window_size, mode='valid')
# idx_strt = int((window_size-1)/2)
# idx_end = int(-((window_size-1)/2))
# d_pt_o = d_pt_o[idx_strt:idx_end]
# d_pt_m2 = d_pt_m2[idx_strt:idx_end]
# d_pt_o_full = d_pt_o_full[idx_strt:idx_end]

# removing names means lines will connect thru gaps
nan_mask = np.isnan(d_o)
d_o = d_o[~nan_mask]
d_pt_o = d_pt_o[~nan_mask]
d_m2 = d_m2[~nan_mask]
d_pt_m2 = d_pt_m2[~nan_mask]
nan_mask = np.isnan(d_o_full)
d_o_full = d_o_full[~nan_mask]
d_pt_o_full = d_pt_o_full[~nan_mask]

# for publication the fill between plots are made using
# 'seasonal' and no moving avg
ax.fill_between(d_pt_o, d_o, d_m2, color='lightgrey', label='', alpha=0.7, zorder=0)
#ax.scatter(dat_pytime, dat_davg, label='Obs', s=1, color='r', marker='o')
#ax.scatter(dat_pytime, mod_davg, label='Mod', s=1, color='g', marker='x')
ax.plot(d_pt_o, d_o, label='Obs.', color='k', linestyle='-', linewidth=1)
# Filter indices of dates before the target date
target_date = datetime(1980, 4, 1)
filtered_indices = [index for index, date in enumerate(d_pt_o_full) if date < target_date]
ax.plot(d_pt_o_full[filtered_indices], d_o_full[filtered_indices], label='', color='k', linestyle='-', linewidth=1)
ax.plot(d_pt_m2, d_m2, label='HOTSSea', color='r', linestyle='--', linewidth=1)
ax.set_ylim(-1.2,1.2)

# Add double-headed arrow and label
ax.axvline(x=datetime(1980, 1, 1), color='k', linestyle='--')
ax.text(datetime(1974, 1, 1), 0.9, 'Cold Period', horizontalalignment='center')
ax.annotate('', xy=(datetime(1968, 1, 1), 0.8),
            xytext=(datetime(1980, 1, 1), 0.8),
             arrowprops=dict(arrowstyle='<->', color='k'))

ax.xaxis.set_major_locator(mdates.YearLocator(5))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.grid(which='major', color='lightgrey')
ax.axhline(y=0, color='black', linestyle='--', linewidth='0.5')

plt.title("HOTSSea v1.02 Modelled vs. Observed Temperature Anomalies at Nanoose Stn, " + str(depth_min) + " to " + str(depth_max) + " m")
plt.ylabel('Anomaly (deg C)')
plt.legend()
plt.savefig('modvobs_nanoose_' + time_inc + '.png', dpi=300)
plt.show()


# seasonal
seas = 2
fig, ax = plt.subplots()
#ax.fill_between(dat_pytime_seas_obs[seas], dat_seas_obs[seas], dat_seas_mod2[seas], color='grey', label='', alpha=0.7, zorder=0)
#ax.scatter(dat_pytime, dat_davg, label='Obs', s=1, color='r', marker='o')
#ax.scatter(dat_pytime, mod_davg, label='Mod', s=1, color='g', marker='x')
ax.plot(d_pt_se_o[seas], d_se_o[seas], label='Nanoose Temp. (4.5 - 400 m)', color='r', linestyle='-', linewidth=3)
# Filter indices of dates before the target date
target_date = datetime(1980, 12, 15)
filtered_indices = [index for index, date in enumerate(d_pt_se_o_full[seas]) if date < target_date]
ax.plot(d_pt_se_o_full[seas][filtered_indices], d_se_o_full[seas][filtered_indices], label='', color='r', linestyle='-', linewidth=3)
#ax.plot(dat_pytime_seas_mod2[seas], dat_seas_mod2[seas], label='Mod', color='g', linestyle='--', linewidth=1)
ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.set_ylim(-1.2,1)
ax.set_xlabel('Time')
ax.set_ylabel('Anomaly (deg C)')


# plot mortalities

# path to salmon survival
import pandas as pd
import xarray as xr
import os
import datetime as dt
import scipy.stats as sp
# mort_p = 'C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/1. Salmon/All Species Survival Exploitation/'
# mort_f = 'chinook_M_EWE_TS Jan 13.csv'
# mort_df = pd.read_csv(os.path.join(mort_p,mort_f), header=[0, 1, 2])
# group_names = mort_df.columns.levels[0][:-1]
mort_p = 'C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/1. Salmon/Chinook Survival Freshwater 2022/MODIFIED'
mort_f = 'salish_M_ocntp_sog_1971-2017.csv'
mort_df = pd.read_csv(os.path.join(mort_p,mort_f), header=[0])

group_names = mort_df.columns[1:]
#mort_dt = mort_df.iloc[:, 0].apply(lambda x: datetime.strptime(str(x), '%Y'))
mort_dt = mort_df.iloc[:, 0]
mort_df_np = np.asarray(mort_dt.values)
mort_pydt = [dt.datetime(x,6,15) for x in mort_df_np]
mort_pydt_np = np.asarray(mort_pydt)

mort_mean = np.nanmean(mort_df.iloc[:,1:], axis=1)
mort_std = np.nanstd(mort_df.iloc[:,1:], axis=1)
# Calculate the z-score based on the confidence level
mort_z = np.abs(np.round(sp.norm.ppf((1 + 0.95) / 2), 2)) # from scipy.stats import norm
# Calculate the upper confidence limit
mort_ucl = mort_mean + mort_z * (mort_std / np.sqrt(len(mort_mean)))
mort_lcl = mort_mean - mort_z * (mort_std / np.sqrt(len(mort_mean)))
# take mean of all mort data
mort_mean = ((mort_mean / 4) - 1)

# since twinx not working, transform the mort data to same -1 to 1 scale
# ax2 = ax.twinx()
# ax.plot(mort_pydt_np, mort_mean, label='SoG Chinook Mort. (Ocean Type)', color='b', linewidth=3)
#
# for grp_name in group_names:
#     mort_dat = mort_df[grp_name].values
#     mort_dat = ((mort_dat / 4) - 1)
#     ax.plot(mort_pydt_np, mort_dat, label='', color='grey', linewidth=1)
# ax.legend()
# ax2.set_ylim(-1,1)
# ax2.set_ylabel('Mort')
# ax2.legend()
# ax2.xaxis.set_major_locator(mdates.YearLocator())
# ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

# plt.title("Season " + str(seas) + ": HOTSS v1.02 Modelled vs. Observed Temperature Anomalies at Nanoose Stn, " + str(depth_min) + " to " + str(depth_max) + " m")
plt.savefig('ChinookMortJustTempsForPPT.png', dpi=300)
plt.show()