import os
import shutil
import metpy.calc as mpcalc
import netCDF4 as nc
import numpy as np
from metpy.units import units


# def getrh(t2m, d2m):
#     #################################################################################################
#     # relative humidity calculation
#     # https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_dewpoint.html
#     # use the 'pint' library bundled with metpy to declare units
#     tair_all_K = t2m.filled() * units.degK
#     d2m_all_K = d2m.filled() * units.degK
#     rh = mpcalc.relative_humidity_from_dewpoint(tair_all_K, d2m_all_K).magnitude
#     rh = np.clip(rh, 0, 1)
#     return rh

def getq(d2m,pressure):
    # specific humidity calculation
    pressurePa = pressure.filled() * units.Pa  # assume pressure at surface is good enough for pressure at 2m
    d2mK = d2m.filled() * units.degK
    #metpy 1.0
    #q2m = mpcalc.specific_humidity_from_dewpoint(pressurePa, d2mK).magnitude
    #metpy 0.12
    q2m = mpcalc.specific_humidity_from_dewpoint(d2mK, pressurePa).magnitude
    return q2m


def get_shifted_data(var, year):
    infile1 = os.path.join(inpath, 'ERA5_SalishSea_9vars_{}.nc'.format(str(year)))
    infile2 = os.path.join(inpath, 'ERA5_SalishSea_9vars_{}.nc'.format(str(year+1)))
    with nc.Dataset(infile1) as ncf:
        tmp = np.zeros(ncf[var].shape) - 999
        tmp[:-1,:,:] = ncf[var][1:,:,:]  # second record becomes first, etc
    if year == year_end:
        # we don't have a 2021 file, duplicate last record as a workaround
        tmp[-1, :, :] = tmp[-2, :, :]
    else:
        with nc.Dataset(infile2) as ncf:
            tmp[-1,:,:] = ncf[var][0,:,:]  # first record from next year is last record of this year
    return tmp


#inpath='/home/mid002/WORK4/SalishSea1500/ECMWF_ERA5/original/'
#outpath='/home/mid002/WORK4/SalishSea1500/ECMWF_ERA5/adjusted_md/'
inpath='C://Users//Greig//Documents//GitHub//NEMO-Salish-Sea-2021//code//'
outpath='C://Users//Greig//Documents//GitHub//NEMO-Salish-Sea-2021//data//forcing/ECMWF//ERA5//'
os.makedirs(outpath,exist_ok=True)
year_start = 1979
year_end = 1985

for year in range(year_start, year_end+1):
    print(year)
    infile = os.path.join(inpath, 'ERA5_SalishSea_9vars_{}.nc'.format(year))
    outfile = os.path.join(outpath, 'ERA5_SalishSea_fixed_{}.nc'.format(year))

    shutil.copy2(infile,outfile)  # copy original file

    # load shifted data
    msr = get_shifted_data('msr', year)
    mtpr = get_shifted_data('mtpr', year)
    msnswrfcs = get_shifted_data('msnswrf', year)  # this should be msnswrf
    msnlwrfcs = get_shifted_data('msnlwrf', year)  # this should be msnlwrf

    with nc.Dataset(outfile,'r+') as ncf:
        # add relative humidity variable
        q = ncf.createVariable('q2m', "float32", ("time", "latitude", "longitude"))
        q.units = "dimensionless"
        q.long_name = "specific_humidity"
        q.missing_value = -999
        q.scale_factor = 1
        q.add_offset = 0

        # compute and store q
        q[:,:,:] = q2m = getq(ncf['d2m'][:, :, :], ncf['sp'][:,:,:])

        # write shifted data
        ncf['msr'][:,:,:] = msr
        ncf['mtpr'][:, :, :] = mtpr
        ncf['msnswrf'][:, :, :] = msnswrfcs
        ncf['msnlwrf'][:, :, :] = msnlwrfcs