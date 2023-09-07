# created by G Oldford Aug 30 2023
# purpose: run trend detection analysis on anomaly time series, cell-wise


# Function to compute the slope of the linear trend
# def compute_slope(x, y):
#    regression = LinearRegression()
#    regression.fit(x, y)
#    return regression.coef_[0]
#
# def compute_slope_GO(x,y):
#    m, c = np.linalg.lstsq(x, y, rcond=None)[0]
#    return m, c
#
# def get_bootstrapped_CI(n_iter, ts, ts_time):
#    values = ts.values[~np.isnan(ts.values)]
#    time = ts_time[~np.isnan(ts.values)]
#    # Fit linear regression to the original data
#    # regression = LinearRegression()
#    # regression.fit(time.reshape(-1, 1), values)
#    A = np.vstack([time, np.ones(len(time))]).T
#    orig_slope, orig_bias = compute_slope_GO(A,values)
#
#    # Initialize an array to store the slopes
#    bootstrap_slopes = np.zeros(n_iter)
#    bootstrap_biases = np.zeros(n_iter)
#
#    # Perform bootstrapping
#    for i in range(n_iter):
#
#        # Generate a bootstrap sample by resampling with replacement
#        indices = np.random.choice(len(values), len(values), replace=True)
#        bootstrap_sample = values[indices]
#        bootstrap_time = time[indices]
#        A = np.vstack([bootstrap_time, np.ones(len(bootstrap_time))]).T
#
#        # Compute the slope of the linear trend on the bootstrap sample
#    #     bootstrap_slope = compute_slope_GO(bootstrap_time.reshape(-1, 1), bootstrap_sample)
#        bootstrap_slope, bootstrap_b = compute_slope_GO(A, bootstrap_sample)
#        bootstrap_slopes[i] = bootstrap_slope
#        bootstrap_biases[i] = bootstrap_b
#
#    return bootstrap_slopes, bootstrap_biases, orig_slope, orig_bias

import os
import time

import numpy as np
import pymannkendall as mk
import xarray as xr

import timeit
from scipy.stats import kendalltau

base_dir = '/project/6006412/goldford/data_temp/extract_results/'
#pattern = 'SalishSea1500-RUN216_1d_grid_T_y{}m{:02d}.nc'

base_dir = 'D:/temp_nemo/trend_output/'
# pattern = 'SalishSea1500-RUN216_1d_grid_T_y{}m{:02d}.nc'

# this is to test whether Mann Kendall gives higher significance
# to weekly sampling rather than seasonal avg anom
use_avg = True # doesn't work - crashes
use_avg_length = "month"
alpha = 0.05


def mann_kendall_trend(data):
    # You may need to preprocess the data here if it's not in the right format
    # For example, convert it to a pandas Series if it's not already
    # data_series = pd.Series(data)

    # Apply the Mann-Kendall test
    result = mk.original_test(data)

    return result.slope, result.h


def mann_kendall_seasonal(data, period):
    # You may need to preprocess the data here if it's not in the right format
    # For example, convert it to a pandas Series if it's not already
    # data_series = pd.Series(data)

    # Apply the Mann-Kendall test
    result = mk.seasonal_test(data, period=period)

    return result.slope, result.h


def trends_3D(variables):
    for variable in variables:
        for depth_group in depth_groups:
            for season in seasons:
                filename = f'weekly_anom_{variable}_{depth_group}_{season}_1980-2018.nc'
                anom_f = os.path.join(base_dir, 'anomalies/', filename)

                print("getting trend for " + variable + " in " + season + " and depth group " + depth_group)
                ds = xr.open_dataset(anom_f)  # Open the NetCDF file using xarray
                print(ds[variable].shape)

                #
                #tau_array = np.apply_along_axis(lambda x: kendalltau(range(len(x)), x)[0], axis=0, arr=ds[variable])
                #print(tau_array)
                # Calculate yearly average for each cell
                #yearly_avgs = ds[variable].groupby("time_counter.year").mean(dim="time_counter")

                # Mann Kendall test
                # pymannkendall
                # takes a vector list as input
                # slopes, trend_tf = xr.apply_ufunc(
                #     mann_kendall_trend,
                #     yearly_avgs,
                #     input_core_dims=[["year"]],
                #     output_core_dims=[[], []],
                #     vectorize=True,
                #     dask="allowed"  # ,
                #     # output_dtypes=[float]
                # )



                ##################
                if use_avg:
                    # Calculate yearly average for each cell
                    yearly_avgs = ds[variable].groupby("time_counter." + use_avg_length).mean(dim="time_counter")

                    # dask setup
                    # chunk_size_time = ds[variable].shape[0] // num_cores
                    # yearly_avgs_dask = da.from_array(yearly_avgs, chunks=(chunk_size_time,))

                    # Mann Kendall test
                    # pymannkendall
                    # takes a vector list as input
                    slopes, trend_tf = xr.apply_ufunc(
                        mann_kendall_trend,
                        yearly_avgs,
                        input_core_dims=[[use_avg_length]],
                        output_core_dims=[[],[]],
                        vectorize=True,
                        dask="allowed"#,
                        #output_dtypes=[float]
                    )

                else:
                    # Mann Kendall test
                    # pymannkendall
                    # takes a vector list as input
                    slopes, trend_tf = xr.apply_ufunc(
                        mann_kendall_trend,
                        ds[variable],
                        input_core_dims=[["time_counter"]],
                        output_core_dims=[[],[]],
                        vectorize=True,
                        dask="allowed"#,
                        #output_dtypes=[float]
                    )

                    # ds_out = results_ds.to_dataset()()
                    # ds_out["TF"] = trend_tf

            # Create a new dataset for writing
            trend_dataset = xr.Dataset(
                {
                    "nav_lat": ds["nav_lat"],
                    "nav_lon": ds["nav_lon"],
                    "slope": (("y", "x"), slopes.values),
                    "sig_tf": (("y", "x"), trend_tf.values),
                },
                coords={"y": ds["y"], "x": ds["x"]},
            )

            # Set attributes for variables
            trend_dataset["nav_lat"].attrs = ds["nav_lat"].attrs
            trend_dataset["nav_lon"].attrs = ds["nav_lon"].attrs

            # Write the dataset to a new NetCDF file
            output_filename = f'anom_trend_{variable}_{depth_group}_{season}.nc'
            print(output_filename)
            if use_avg:
                output_path = os.path.join(base_dir, 'anom_trend_frommonthly/', output_filename)
            else:
                output_path = os.path.join(base_dir, 'anom_trend_fromweekly/', output_filename)
            trend_dataset.to_netcdf(output_path)

            ds.close()  # Close the dataset


def trends_2D(variables):
    for variable in variables:
        for season in seasons:
            filename = f'weekly_anom_{variable}_{season}_1980-2018.nc'
            anom_f = os.path.join(base_dir, 'anomalies/', filename)

            print("getting trend for " + variable + " in " + season)
            ds = xr.open_dataset(anom_f)  # Open the NetCDF file using xarray

            # Calculate yearly average for each cell
            yearly_avgs = ds[variable].groupby("time_counter.year").mean(dim="time_counter")

            #trend through each pixel-slice
            # https://stackoverflow.com/questions/66594056/linear-regression-on-each-grid-cell-across-time-dim
            lin_fit = yearly_avgs.polyfit('year', deg=1, skipna=True)
            # y = mx + b
            a = lin_fit.sel(degree=1) # slopes of linregress
            b = lin_fit.sel(degree=0)

            # Create a new dataset for writing
            trend_dataset = xr.Dataset(
                {
                    "nav_lat": ds["nav_lat"],
                    "nav_lon": ds["nav_lon"],
                    "slope": (("y", "x"), a['polyfit_coefficients'].values),
                },
                coords={"y": ds["y"], "x": ds["x"]},
            )

            # Set attributes for variables
            trend_dataset["nav_lat"].attrs = ds["nav_lat"].attrs
            trend_dataset["nav_lon"].attrs = ds["nav_lon"].attrs

            # Write the dataset to a new NetCDF file
            output_filename = f'trend_slopes_{variable}_{season}.nc'
            output_path = os.path.join(base_dir, 'trend_output/', output_filename)
            trend_dataset.to_netcdf(output_path)

            ds.close()  # Close the dataset


years = np.arange(1980, 2019, 1)
seasons = ['winter', 'spring', 'summer', 'fall']
depth_groups = ["0to30m", "30to150m", "gt150m", "allz"]
variables = ['votemper', 'vosaline', 'vomecrty', 'vozocrtx']
seasons = ['summer', 'fall']
variables = ['votemper']
depth_groups = ["allz"]

if __name__ == "__main__":

    print(time.localtime())
    trends_3D(variables)
    print(time.localtime())
    #print(execution_time)
    #variables = ['mldkz5', 'mldr10_1']
    # trends_2D(variables)
