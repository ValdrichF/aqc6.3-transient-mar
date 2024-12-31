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

ibound = ibound.sel(time='20110926', drop=True)
x = ibound.coords["x"]
y = ibound.coords["y"]
layer = ibound.coords["layer"]

# import the original recharge file to add to the newly created one
# org_rch = imod.idf.open(model_path / "GRONDWATERAANVULLING//gemiddelde_grondwateraanvulling_new_25.IDF")

sample = pd.read_csv(r"D:\Dropbox\WUR\Rscripts\paper_2\ortho_grid.csv")
grp_no = sample.scenario.unique()
rch_no = sample.rch_site.unique()

save_dir = (model_path / "rch_transient_individual_sites")
(save_dir/ "nc").mkdir(exist_ok=True, parents=True)

# Loop thorugh each group
for grp in grp_no:
	sample_grp = sample.loc[sample['scenario']==grp] # _no
	# Create a dummy rch file with 0s that will be filled in later
    ## Created like one of the idf files
	rch_tmp = xr.zeros_like(ibound)
	# Then through each rch_no
	for row in sample_grp.itertuples(): #_rch
		selection =(
			(layer == 1)
			& (y   >= row.S_side)
			& (y   <= row.N_side)
			& (x   >= row.W_side)
 			& (x   <= row.E_side) # sample_rch.iloc[0]
			)
		rch_tmp = xr.where(selection, row.recharge, rch_tmp) # otherwise use sample_rch.iloc[0].recharge
		
	imod.idf.save(save_dir / f"rch_{grp}.idf", rch_tmp) # _{rech} # +org_rch
	rch_tmp = rch_tmp.sel(layer = 1,drop = True)
	rch_tmp = rch_tmp.to_dataset(name = "rch")
	rch_tmp = rch_tmp.astype(np.float32)
	rch_tmp.to_netcdf(save_dir / "nc" / f"rch_{grp}.nc")

        
        
        
        
        
        
        
