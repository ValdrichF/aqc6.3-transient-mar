# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 18:24:29 2023

@author: valdf

Extract the response for each cell to a CSV
Filter it to record only the response > 1cm
Also record the scenario and rch that it is related to
"""

import imod
import xarray as xr
from pathlib import Path
import pandas as pd

res_dir = Path(r"Z:\MODFLOW_transient_res")
big_grid = pd.read_csv(r"C:\Users\Valdrich\Dropbox\WUR\Rscripts\paper_2\ortho_grid.csv") 

scenarios = [n for n in res_dir.iterdir() if n.is_dir()]
del scenarios[0:3]

# Writing to the csv within each loop
# delete the csv if it exists
if (res_dir/"response_cells_thresh.csv").exists():
	(res_dir/"response_cells_thresh.csv").unlink()

header = True

for flux_dir in scenarios:
	print(f"scenario: {flux_dir.name}")
	
	flux_dir = flux_dir/'response'
	
	idf = imod.idf.open(flux_dir/"*L1.idf")
	
	# Only where there's a response >1cm
	idf = idf.where(idf>=0.01)
	scenario_dir = list(flux_dir.parents)[-3]
	
	nearest_rch = imod.idf.open(scenario_dir/"nearest_rch"/"rch_l1.idf") # Recharge rate
	
	idf = idf.sel(layer=1, drop=True)
	nearest_rch = nearest_rch.sel(layer=1, drop=True)
	
	idf_interest = idf.where(nearest_rch.notnull())

	scenario_num = int(scenario_dir.name.split("_")[1])
	
	LHS_subset = big_grid.loc[(big_grid["scenario"] == scenario_num)]
	
	rch_site = xr.zeros_like(nearest_rch)
	
	for row in LHS_subset.itertuples():
		
		rch_site = xr.where(nearest_rch==row.recharge, row.rch_site, rch_site)
		
	ds = idf_interest.to_dataset(name = "response")
	ds["rch_site"] = rch_site
	ds["scenario"] = scenario_num
	
	table = ds.to_dataframe()
	table = table[table.notnull().all(1)]
	table.to_csv(res_dir/"response_cells_thresh.csv", mode = 'a', header = header)
	
	header = False
	
