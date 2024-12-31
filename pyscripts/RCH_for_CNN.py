# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 15:15:54 2022

@author: Valdrich

Creating the recharge files:
loop through the csv files
loop through each row in the csv and make the recharge file

Run this file after running the Scenarios_designOALHS.R file
"""
# import the packages
import imod 
import xarray as xr
import pandas as pd
import numpy as np

from pathlib import Path

model_path = r"F:\Modeldatabase\AMIGO_needed"
model_path = Path(model_path)
ibound = imod.idf.open(model_path / "SHD"/"head_20110926_l1.idf")

ibound = ibound.sel(time='20110926', layer = 1, drop=True)
x = ibound.coords["x"]
y = ibound.coords["y"]

#Crop the model to the catchment rectangle
# Also need to pad 0s to help design the ML model
h = 832
w = 1472

# Find the number of cells in the x and y directions
Xmin = 204519
Xmax = 241200

Ymin = 439298
Ymax = 460020

Ny = ibound.sel(x = slice(Xmin, Xmax), y = slice(Ymax, Ymin)).shape[0]
Nx = ibound.sel(x = slice(Xmin, Xmax), y = slice(Ymax, Ymin)).shape[1]
dx = ibound.dx
dy = ibound.dy

top_pad   = np.floor((h - Ny) / 2).astype(np.uint16)
bottom_pad = np.ceil((h - Ny) / 2).astype(np.uint16)
right_pad = np.ceil((w - Nx) / 2).astype(np.uint16)
left_pad = np.floor((w - Nx) / 2).astype(np.uint16)

ibound = ibound.sel(x = slice(Xmin-left_pad*dx, Xmax+right_pad*dx),
                    y = slice(Ymax-top_pad*dy, Ymin+bottom_pad*dy))

ibound[:top_pad] = 0
ibound[slice(-1*bottom_pad,None)] = 0

ibound[:,:left_pad] = 0
ibound[:,slice(-1*right_pad,None)] = 0
x = ibound.coords["x"]
y = ibound.coords["y"]


# import the original recharge file to add to the newly created one
# org_rch = imod.idf.open(model_path / "GRONDWATERAANVULLING//gemiddelde_grondwateraanvulling_new_25.IDF")

sample = pd.read_csv(r"D:\Dropbox\WUR\Rscripts\paper_2\ortho_grid.csv")
grp_no = sample.scenario.unique()
rch_no = sample.rch_site.unique()

save_dir = (model_path / "rch_transient_individual_sites")
(save_dir/ "nc").mkdir(exist_ok=True, parents=True)

# Loop thorugh each group
# Create a dummy rch file with 0s that will be filled in later
## Created like one of the idf files
rch = xr.zeros_like(ibound)
# Then through each rch_no
for row in sample.itertuples(): #_rch
	selection =(
		(y   >= row.S_side)
		& (y   <= row.N_side)
		& (x   >= row.W_side)
		& (x   <= row.E_side) # sample_rch.iloc[0]
		)
	rch_tmp = xr.where(selection, row.recharge, rch) # otherwise use sample_rch.iloc[0].recharge
	
	f_name = f"rch_{row.scenario}_{row.rch_site}"
	print(f"Working on {f_name}")
	imod.idf.save(save_dir / (f_name + ".idf"), rch_tmp) # _{rech} # +org_rch
	# rch_tmp = rch_tmp.sel(layer = 1,drop = True)
	rch_tmp = rch_tmp.to_dataset(name = "rch")
	rch_tmp = rch_tmp.astype(np.float32)
	rch_tmp.to_netcdf(save_dir / "nc" / (f_name + ".nc"))

        
        
        
        
        
        
        
