import netCDF4 as nc
import os
import re

# Necessary defaults
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import os, sys

# Xarray
import xarray as xr
# Dask stuff
import dask.array as da
from dask.diagnostics import ProgressBar

from numba import jit  # Speedup for python functions

# ----------------------
# Parameters
dates = pd.date_range('1984-01-01', periods=34, freq="A").year


@jit(nogil=True)
def mk_cor(x, y, pthres=0.05, direction=True):
    # Check NA values
    co = np.count_nonzero(~np.isnan(x))
    if co < 4:
        return -9999
    # The test
    tau, p_value = stats.kendalltau(x, y)

    # Criterium to return results
    if p_value < pthres:
        # Check direction
        if direction:
            if tau < 0:
                return -1
            elif tau > 0:
                return 1
        else:
            return tau
    else:
        return 0


def kendall_correlation(x, y, dim='year'):
    return xr.apply_ufunc(
        mk_cor, x, y,
        input_core_dims=[[dim], [dim]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[int]
    )


def save_xarray(outpath, arr):
    arr.to_netcdf(os.path.join(outpath,),
                  encoding={'trend': {'dtype': 'int16',
                                      'zlib': True,
                                      'complevel': 9,
                                      '_FillValue': -9999}})


def plot_route(route, id, year):
    route.plot(robust=True, cmap='viridis')
    plt.title('The route %s in year %s' % (id, year))
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')


def process_route(route, outpath):


    # Open the NetCDF file
    with nc.Dataset(route) as dataset:
        # Access and manipulate the data as needed
        # For example, you can retrieve variables using dataset.variables['variable_name']
        # and perform operations on the data
        ds = xr.Dataset({'EVI': (['x', 'y', 'year'], dataset['EVI'][:])})
        ds.coords["year"] = dates
        ds = ds.rename({'x': 'Longitude', 'y': 'Latitude', 'year': 'year'})

        print("Starting to process %s GB of data" % (round(ds.nbytes * (2 ** -30), 2)))
        ds2 = ds.pipe(lambda x: x * 0.0001).where(ds.EVI > 0)
        x = xr.DataArray(np.arange(len(ds2['year'])) + 1, dims='year', coords={'year': ds2['year']})
        with ProgressBar():
            r = kendall_correlation(ds2, x, 'year').compute()
        r = r.rename({'EVI': 'trend'})
        save_xarray(outpath, r)


# Function
if __name__ == '__main__':
    print("Start processing")

    # Testing
    route = "D:/temp_nemo/trend_output/anomalies/weekly_anom_votemper_30to150m_spring_1980-2018.nc"
    outpath = "D:/temp_nemo/trend_output/temp/anom_trend_votemper_0to30m_winter.nc"

    process_route(route, outpath)