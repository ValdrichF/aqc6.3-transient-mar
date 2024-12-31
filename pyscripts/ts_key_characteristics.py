# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:16:22 2023

@author: Valdrich

output - csv of total, maximum and area of the response


"""
from pathlib import Path
import imod
import pandas as pd
from matplotlib import pyplot as plt
import xarray as xr


res_dir = Path(r"Z:\MODFLOW_transient_res")
big_grid = pd.read_csv(r"D:\Dropbox\WUR\Rscripts\paper_2\ortho_grid.csv")
sample = imod.idf.open(res_dir/"rch_1_l1/response/*L1.idf")

sy = imod.idf.open(res_dir/"avg_baseline/metaswap/msw_sc1/msw*l1.idf") 
sy = sy.sel(layer=1, drop=True)
sy = sy.reindex_like(sample.time, method = 'nearest')
sy = sy.load()


def create_ts_flux(flux_dir):
	print(f"scenario: {flux_dir.parent.name}")
	
	idf = imod.idf.open(flux_dir/"*L1.idf")
	idf = idf.sel(layer=1, drop=True)
	# Only where there's a response >1cm
	idf = idf.where(idf>=0.01, 0)
	vol_water = idf*sy
	ds = xr.Dataset({
		"total":idf,
		"total_water":vol_water
		})

	scenario_dir = list(flux_dir.parents)[-3]
	scenario_num = int(scenario_dir.name.split("_")[1])
	LHS_subset = big_grid.loc[(big_grid["scenario"] == scenario_num)]
	
	# split into individual sites
	nearest_rch = imod.idf.open(scenario_dir/"nearest_rch"/"rch_l1.idf")
	nearest_rch = nearest_rch.sel(layer=1, drop=True)
	nearest_rch3D = xr.zeros_like(idf, dtype = bool).expand_dims(rch_site = LHS_subset.rch_site)
	for rch_site in  LHS_subset.rch_site: # _no
		nearest_rch3D.loc[rch_site,:,:] = nearest_rch==rch_site
	ds_individual = ds.where(nearest_rch3D)


	total = ds_individual.sum(dim = ["x", "y"])*ds_individual.dx**2
	total["maximum"] = ds_individual.total.max(dim = ["x", "y"])
	total["area"] = (ds_individual.total>0).sum(dim = ["x", "y"])*ds_individual.dx**2
	
	# asign scenario number as coordinate
	total = total.assign_coords({"scenario":scenario_num}).load()
	return total


scenarios = [n for n in res_dir.iterdir() if n.is_dir()]
del scenarios[0:3]

df = [create_ts_flux(n/'response') for n in scenarios]
df = xr.concat(df, dim = "scenario")
df = df.to_dataframe()
df.to_csv(res_dir/"response_sites_thresh0.csv")
df.to_csv(r"D:\Dropbox\WUR\Rscripts\paper_2\response_sites_thresh0.csv")

if False:
	a = sy.isel(
		time=slice(0, 53, 4)
		).plot.pcolormesh(
			levels = [0, 0.1, 0.2, 0.3], col = "time", col_wrap=4, aspect = sy.x.shape[0]/sy.y.shape[0]
			)
	axes = a.axs
	
	for row in axes:
		for ax in row:
			ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
	plt.savefig(res_dir / "storage_coefficient_baseline.png", dpi = 300)
	plt.draw()
	plt.clf()

	a = sy.mean(
		dim = 'time'
		).plot.pcolormesh(
			levels = [0, 0.1, 0.2, 0.3], aspect = sy.x.shape[0]/sy.y.shape[0], size = 3
			)
	axes = a.axes
	axes.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
	plt.savefig(res_dir / "mean_storage_coefficient_baseline.png", dpi = 300)
	plt.draw()
	plt.clf()

