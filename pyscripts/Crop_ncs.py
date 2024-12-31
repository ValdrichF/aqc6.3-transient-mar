# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 12:37:04 2024

@author: valdf
"""

# crop the nc files
import xarray as xr
from pathlib import Path

nc_dir = Path(r"F:\Modeldatabase\AMIGO_needed\rch_transient\nc")
ncs = list(nc_dir.glob("*.nc"))
save_dir = nc_dir.parent / "cropped"
save_dir.mkdir(exist_ok=True)

common_path = Path(r"D:\Modeldatabase\common_ncs_New")
common = xr.open_dataset(common_path / "common_data_DE_KD_VR_DC_DL.nc")

for n in ncs:
	
	# read the nc
	nc_file = xr.open_dataset(n)
	
	# Crop to the extent
	nc_file_cropped = nc_file.reindex_like(common)
	
	# save
	nc_file_cropped.to_netcdf(save_dir / n.name)








