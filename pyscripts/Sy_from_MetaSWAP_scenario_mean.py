# -*- coding: utf-8 -*-
"""
Created on Wed May  8 11:50:23 2024

@author: valdf

dh weighted mean storage yield from the MetaSWAP
	- find absolute dh for each stress period
	- the storage yield from MetaSWAP for that section of dh
	- multiply the dh with the specific yield divided by the sum dh
"""

from pathlib import Path
import xarray as xr
import imod
import numpy as np

sy_scenario = Path(r"Z:\Sy_scenario\results")

heads = imod.idf.open(sy_scenario/"head/head*.idf")
heads = heads.sel(layer = 1)
SHD = imod.idf.open(sy_scenario.parent/"SHD/SHD*.idf")
SHD = SHD.sel(layer = 1)
# add time dimension for SHD
init_time = heads.time[0].values-np.timedelta64(7, 'D')
SHD = SHD.assign_coords(time = init_time)
SHD = SHD.expand_dims(dim="time")

# merge SHD with the heads
heads = xr.concat([SHD, heads], dim = "time")

# specific yield:
sy = imod.idf.open(sy_scenario/"metaswap/msw_S01/msw_S01*.IDF")
sy = sy.sel(layer = 1)
sy = sy.reindex_like(heads.time.isel(time = slice(1, None)),
					 method = "nearest")

########### diff head weighted #####
dhead = heads.diff(dim = "time")
dhead = dhead.where(dhead>0, 0)

weighted_sy = (sy*dhead).sum(dim="time")/dhead.sum(dim="time")
imod.idf.save(sy_scenario.parent/"mean_sy", weighted_sy)




