# -*- coding: utf-8 -*-
"""
Created on Mon May  6 11:09:59 2024

@author: valdf

MetaSWAP-iMODFLOW scenario
- low initial head
- constant recharge of 10 mm/day
- constant evaporation of .1 mm/day (try to skip)
- scenario duration of 90 days (Done in R)
"""

import imod
from pathlib import Path
import xarray as xr

baseline_heads = imod.idf.open(r"Z:\MODFLOW_transient_res\avg_complete\head\*.idf")
baseline_heads = baseline_heads.isel(time = slice(53,None))

SHD = baseline_heads.min(dim = 'time')

imod.idf.save(r"Z:\Sy_scenario\SHD\SHD", SHD)

# Only recorded the heads for Layers 1, 4, 5, 7
# The other layers have similar heads as it, so copy it and rename them appropriately
# The scenario 'average_initial' has the heads for all the layers.

# Recharge file
RCH = xr.ones_like(SHD.sel(layer = 1))
RCH = RCH*10
imod.idf.save(r"Z:\Sy_scenario\RCH\RCH", RCH)

# Recharge file
ET = xr.zeros_like(SHD.sel(layer = 1))
imod.idf.save(r"Z:\Sy_scenario\ET\ET", ET)



