# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 12:31:48 2024

@author: valdf

Compare the stationary response to the transient response at the end of recharge (25-02-2012)

"""
import imod
from pathlib import Path
import contextily as cx
from matplotlib import pyplot as plt
import xarray as xr
import pandas as pd
import re
from functions import dist_to_rch
from joblib import Parallel, delayed


# Combine the stationary responses for each group
big_grid = pd.read_csv(r"D:\Dropbox\WUR\Rscripts\paper_2\ortho_grid.csv")
# big_grid = pd.read_csv(r"D:\Dropbox\WUR\Rscripts\test_ml_grid_10ha.csv")
# preds_dir = Path(r"Z:\Stationary_preds_multiple\preds_UNET_1000\Predictions")

# Adding river data
p = Path(r"D:\Modeldatabase\common_ncs_New\common_data_DE_KD_VR_DC_DL.nc")
ds = xr.open_dataset(p)
riv_stage = ds.depth-ds.level
riv_conductance = ds.conductance
riv_conductance = riv_conductance.where(riv_conductance>0.35)
riv_conductance = riv_conductance.where(riv_stage<=1)
riv_conductance = riv_conductance.where(riv_stage!=0)

r = xr.ones_like(ds.depth).where(riv_conductance.notnull())



#%% Check!!!!
path_A = Path(r"Z:\MODFLOW_transient_res")
path_B = Path(r"Z:\Stationary_res_AMIGO\idfs")

results_dir = path_B.parent
(results_dir/"png_divergent").mkdir(exist_ok=True, parents=True)

groups = path_B.glob("*.idf")

# Specific yield
sample = imod.idf.open(path_A/"rch_1_l1/response/*L1.idf")
sy = imod.idf.open(path_A/"avg_baseline/metaswap/msw_sc1/msw*l1.idf") 
sy = sy.sel(layer=1, drop=True)
sy = sy.reindex_like(sample.time, method = 'nearest')
sy = sy.load()

# Legend
# standard colors for sequential and divergent legends
sequential_col = ['#feedde', '#fdd0a2', '#fdae6b', '#fd8d3c', '#f16913', '#d94801', '#8c2d04']
divergent_col =  ['#e66100', '#eb8a54', '#e7b195', '#D3D3D3', '#b0a1c4', '#886cb0', '#5d3a9b']
label_1 = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32]
label_2 = [-1, -0.5, -0.05, 0.05, 0.5, 1]
label_3 = [-0.1, -0.05, -0.01, 0.01, 0.05, 0.1]
label_4 = [0.01, 0.2, 0.4, 0.6, 0.8, 0.99] 

AB_label = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64]
AB_col = ['#ffffff','#dbe1e8','#b8c4d2','#94a6bb','#7089a5','#4c6b8e','#294e78', '#053061']


C_label = label_3
C_col = divergent_col

D_label = label_2
D_col = divergent_col

#%%

# Some stuff for the figures
source = cx.providers.CartoDB.Voyager
kwargs_basemap = {"zoom":10, "attribution_size":2, "attribution":""}

# m to cm
AB_label = [n*100 for n in AB_label]
C_label = [n*100 for n in C_label]
D_label = [n*100 for n in D_label]

# reverse the colors
AB_col = AB_col[::-1]
C_col = C_col[::-1]
D_col = D_col[::-1]

# labels for the figures so they can be referenced as figure 5A
labels = [["A", "B"], ["C", "D"]]


totals = []
max_list = []
area_list = []
for group in groups:
	
	group = group.stem
	print(f"working on files starting with {group}")
	
	grp_no = re.search(r'\d+', group).group()
	grp_no = int(grp_no)

	#%% CHECK!!!
	sample_grid = big_grid.loc[(big_grid['scenario']==grp_no)]
	preds = imod.idf.open(path_B / (group + ".idf")) # Check!!!
	scenario = imod.idf.open(path_A / group / "response/response_20120225000000_l1.idf") # Check!!!
	scenario = scenario.isel(time = 0, layer = 0,drop = True)
	preds = preds.sel(layer = 1,drop = True)
	
	#%%%
	
	# Set minimum to reduce the maximum ratio
	scenario = scenario.where(scenario>0.01, 1e-8).fillna(1e-8)
	preds = preds.where(preds>0.01, 1e-8)

	difference = preds-scenario
	ratio = (preds-scenario)*2/(preds+scenario)

	# Set minimum to reduce the maximum ratio
	scenario = scenario.where(scenario>0.01, 1e-8).fillna(1e-8)
	preds = preds.where(preds>0.01, 1e-8)
	
	# Calculate the total response for each recharge site
	nearest_rch = imod.idf.open(path_A/group/"nearest_rch"/"rch_l1.idf")
	nearest_rch = nearest_rch.isel(layer=0, drop=True)
	
	# split into individual sites
	nearest_rch3D = xr.zeros_like(
		scenario, dtype = bool
		).expand_dims(rch_site = sample_grid.rch_site)
	for rch_site in  sample_grid.rch_site: # _no
		nearest_rch3D.loc[rch_site,:,:] = nearest_rch==rch_site

	# split the combined dataset (preds+scenario) by the recharge sites
	ds = xr.Dataset({"Transient":scenario,
				    "Stationary":preds})
	ds_individual = ds.where(nearest_rch3D)
	
	# total
	ds_total = ds_individual.sum(dim=['x', 'y'])*ds_individual.dx**2
	ds_total = ds_total.assign_coords({"scenario":grp_no}).load()
	totals.append(ds_total.load())
	# maximum
	ds_max = ds_individual.max(dim=['x', 'y'])
	ds_max = ds_max.assign_coords({"scenario":grp_no}).load()
	max_list.append(ds_max.load())
	# area
	ds_area = (ds_individual>0.01).sum(dim=['x', 'y'])*ds_individual.dx**2
	ds_area = ds_area.assign_coords({"scenario":grp_no}).load()
	area_list.append(ds_area.load())
	
totals = xr.concat(totals, dim = 'scenario')
totals_csv = totals.to_dataframe()
totals_csv.to_csv(results_dir/"total_response.csv")

max_list = xr.concat(max_list, dim = 'scenario')
max_csv = max_list.to_dataframe()
max_csv.to_csv(results_dir/"max_response.csv")

area_list = xr.concat(area_list, dim = 'scenario')
area_csv = area_list.to_dataframe()
area_csv.to_csv(results_dir/"area_response.csv")


def plot_model_diff(group):

	group = group.stem
	print(f"working on files starting with {group}")
	
	#%% CHECK!!!
	preds = imod.idf.open(path_B / (group + ".idf")) # Check!!!
	scenario = imod.idf.open(path_A / group /("response/response_20120225000000_l1.idf")) # Check!!!
# 	scenario = imod.idf.open(path_A / (group + ".idf")) # Check!!!
	scenario = scenario.isel(layer = 0, time = 0,drop = True)
	preds = preds.sel(layer = 1,drop = True)
	#%%
	
	# Set minimum to reduce the maximum ratio
	scenario = scenario.where(scenario>0.02, 1e-8).fillna(1e-8)
	preds = preds.where(preds>0.02, 1e-8)

	difference = preds-scenario
	ratio = (preds-scenario)/(preds+scenario)*2

	# Transparent background
	difference = difference.where((preds+scenario)>0.01)
	ratio = ratio.where((preds+scenario)>0.01)
	scenario = scenario.where(scenario>0.01)
	preds = preds.where(preds>0.01)
	
	if True:
		
		
		fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, figsize = (8,4.5))
		imod.visualize.spatial.plot_map(scenario*100, AB_col, AB_label, basemap=source, fig = fig, ax = axes[0,0], kwargs_basemap=kwargs_basemap)
		imod.visualize.spatial.plot_map(preds*100, AB_col, AB_label, basemap=source, fig = fig, ax = axes[0,1], kwargs_basemap=kwargs_basemap)
		imod.visualize.spatial.plot_map(difference*100, C_col, C_label, basemap=source, fig = fig, ax = axes[1,0], kwargs_basemap=kwargs_basemap)
		imod.visualize.spatial.plot_map(ratio*100, D_col, D_label, basemap=source, fig = fig, ax = axes[1,1], kwargs_basemap=kwargs_basemap)
		
		for row, label_row in zip(axes, labels):
			for ax, label in zip(row, label_row):
				ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
				ax.set_title(label, loc = 'left', 
							 fontdict = {
								 'fontsize':12,
								 'fontweight':'bold'								 
								 })
				r.plot(ax = ax, add_colorbar=False, add_labels = False, cmap = 'PuBuGn')

		plt.tight_layout()
		plt.savefig(results_dir / f"png/rch{group[3:]}.png", dpi = 300)
		plt.draw()
		plt.clf()
		
# 	res_ds = xr.Dataset({
# 		'Transient':scenario,
# 		'Steady':preds,
# 		'difference':difference,
# 		'ratio':ratio
# 		})
 	
# 	response_rivers = xr.merge([res_ds, ds])
# 	response_rivers = response_rivers.where(riv_stage!=0)
# 	response_rivers = response_rivers.assign_coords(group = group)
# 	return response_rivers


groups = path_B.glob("*.idf")
responses = [plot_model_diff(n) for n in groups]
responses = xr.concat(responses, dim = 'group')
responses_df = responses.to_dataframe()
responses_df = responses_df.dropna()
responses_df.to_csv(results_dir/"response_rivers.csv")

# Parallel(n_jobs=24)(delayed(plot_model_diff)(n) for n in groups)







