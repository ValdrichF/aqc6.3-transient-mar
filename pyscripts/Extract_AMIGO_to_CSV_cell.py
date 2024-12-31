# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 14:18:09 2023

@author: Valdrich

output - csv of avg KD, VCW, and the depth to groundwater within 1km from the cell

	Drain and Rivers, only consider those that are active, 
	ie stage below the groundwater level or infiltration factor == 1
	
"""
from pathlib import Path
import imod
import pandas as pd
import xarray as xr
from scipy.ndimage import gaussian_filter
from dask import delayed, compute
import numpy as np

res_dir = Path(r"Z:\MODFLOW_transient_res")
big_grid = pd.read_csv(r"D:\Dropbox\WUR\Rscripts\paper_2\ortho_grid.csv") 
# big_grid = pd.read_csv(r"C:\Users\Valdrich\Dropbox\WUR\Rscripts\paper_2\ortho_grid.csv") 

# maximum baseline head
baseline_heads = imod.idf.open(res_dir/"avg_baseline/head/*l1.idf")
baseline_heads = baseline_heads.sel(layer = 1, drop=True)

# GHG
start = baseline_heads.time.min().values
end = baseline_heads.time.max().values
dates = pd.date_range(start=start, end=end, freq="SMS") + pd.DateOffset(days=13)
head_bimonthly = baseline_heads.reindex(time=dates, method="nearest", tolerance=pd.Timedelta(days=7))
heads_mean = head_bimonthly.rolling(time = 3, center = True).mean()
GHG = heads_mean.max(dim = 'time')
ind_max = abs(baseline_heads-GHG).mean(dim = ["x", "y"]).argmin().values

# KD, VCW
small_model = Path(r"D:\Modeldatabase\Baaksebeek\Small_model")
KD = imod.idf.open(small_model/"Transmissivity/Transmissivity_l1.idf")
VCW = imod.idf.open(small_model/"Resistance/Resistance_l1.idf")
KD = KD.sel(layer = 1, drop = True)
VCW = VCW.sel(layer = 1, drop = True)

# Depths
depths = imod.idf.open(res_dir/"avg_baseline/depth/*l1.idf")
depths = depths.sel(layer = 1, drop=True)

# Storage coefficient of the entire unsaturated zone: from MetaSWAP scenario
sy_scenario = Path(r"Z:\Sy_scenario")
sy = imod.idf.open(sy_scenario/"mean_sy*l1.idf")
sy = sy.sel(layer = 1,drop = True)
# fill nans (rivers or open water) with 1
sy = sy.where(sy.notnull(), 1)
# Saturated water content: from the inputs to MetaSWAP
por = imod.idf.open(sy_scenario/"saturated_water_content.idf")
por = por.where(por.notnull(), por.mean().values)
# Storage coefficient from the baseline scenario
sy_baseline = imod.idf.open(res_dir/"avg_baseline/metaswap/msw_sc1/msw*l1.idf")
sy_baseline = sy_baseline.sel(layer = 1, drop = True)

# Active drains and river
model_dat = Path(r"D:\Modeldatabase\AMIGO31")

legger_stage = imod.idf.open(model_dat/"OPPERVLAKTEWATER/LEGGER/peil_legger_winter.idf")
IJssel_stage = imod.idf.open(model_dat/"OPPERVLAKTEWATER/RIVIER/maand/*.idf")
IJssel_stage = IJssel_stage.sel(time = slice("20111001", "20121001"))
RIV_stage = [imod.idf.open(n) for n in (model_dat/"OPPERVLAKTEWATER/TOP10").glob("Peil*.idf")]
RIV_stage.extend([legger_stage, IJssel_stage])
RIV_stage = xr.concat(RIV_stage, pd.Index([1]*len(RIV_stage), name = "layer"))
RIV_stage = RIV_stage.min(dim="layer")
RIV_stage = RIV_stage.sel(time = baseline_heads.time, method = "ffill")
RIV_stage["time"] = baseline_heads.time
RIV_stage = RIV_stage-baseline_heads

legger_cond = imod.idf.open(model_dat/"OPPERVLAKTEWATER/LEGGER/conductance_legger_winter.idf")
IJssel_cond = imod.idf.open(model_dat/"OPPERVLAKTEWATER/RIVIER/Conductance_ISG_AMIGO.idf")
RIV_cond = [imod.idf.open(n) for n in (model_dat/"OPPERVLAKTEWATER/TOP10").glob("conductance*.idf")]
RIV_cond.extend([legger_cond, IJssel_cond])
RIV_cond = xr.concat(RIV_cond, pd.Index([1]*len(RIV_cond), name = "layer"))
RIV_cond = RIV_cond.max(dim="layer")

drainage = [imod.idf.open(n) for n in (model_dat/"DRAINAGE").glob("peil*.idf")]
ontgrondingen = imod.idf.open(model_dat/"OPPERVLAKTEWATER/ONTGRONDINGEN/peil_ONTGRONDINGEN.idf")
steengroeve = imod.idf.open(model_dat/"DRAINAGE/drainage_niveau_steengroeve.IDF")
drainage.extend([ontgrondingen, steengroeve])
drainage = xr.concat(drainage, pd.Index([1]*len(drainage), name = "layer"))
drainage = drainage.min(dim="layer")
drainage = drainage-baseline_heads
# imod.idf.save(model_dat/"min_drain_level.idf", drainage)

drainage_cond = [imod.idf.open(n) for n in (model_dat/"DRAINAGE").glob("conductance*.idf")]
ontgrondingen_cond = imod.idf.open(model_dat/"OPPERVLAKTEWATER/ONTGRONDINGEN/conductance_ontgrondingen.idf")
steengroeve_cond = imod.idf.open(model_dat/"DRAINAGE/drainage_conductance_steengroeve.IDF")
drainage_cond.extend([ontgrondingen_cond, steengroeve_cond])
drainage_cond = xr.concat(drainage_cond, pd.Index([1]*len(drainage_cond), name = "layer"))
drainage_cond = drainage_cond.min(dim="layer")

# approximate river density, normally drain density, but AMIGO
riv_density = RIV_stage.notnull().astype('float32')

# Select the cells with a river in that stress-period
RIV_cond = RIV_cond.where(riv_density)

# Density of the rivers active during high groundwater
riv_flux = imod.idf.open(res_dir / "avg_baseline/bdgriv/*idf").sel(layer = 1, drop = True)
riv_flux = xr.where(riv_flux.isnull(), 0, riv_flux)

# Where is it non-zero?
riv_flux = np.sign(riv_flux)
riv_flux = abs(riv_flux)
riv_flux = riv_flux.isel(time = ind_max, drop = True)

# Filling nan conductace with 0s
# RIV_cond = RIV_cond.fillna(0)
# drainage_cond = drainage_cond.fillna(0)
# # Fill nan rivers and drains with depth
# RIV_stage = RIV_stage.fillna(depths)
# drainage = drainage.fillna(depths)

# Crop and reshape the data
Xmin = KD.x.min().values
Xmax = KD.x.max().values
Ymin = KD.y.min().values
Ymax = KD.y.max().values

RIV_stage = RIV_stage.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
RIV_cond = RIV_cond.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
drainage = drainage.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
drainage_cond = drainage_cond.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
depths = depths.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
riv_density = riv_density.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
riv_flux = riv_flux.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
sy = sy.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
sy_baseline = sy_baseline.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")
por = por.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin)).transpose(..., "y", "x")

# selecting the winter/spring period for the initial response and the decay rate
depth_winter = depths.sel(time = slice(None, "20120301")).mean(dim = 'time')
RIV_stage_winter = RIV_stage.sel(time = slice(None, "20120301")).mean(dim = 'time')
drainage_winter = drainage.sel(time = slice(None, "20120301")).mean(dim = 'time')
riv_density_winter = riv_density.sel(time = slice(None, "20120301")).mean(dim = 'time')
sy_baseline_winter = sy_baseline.sel(time = slice(None, "20120301")).mean(dim = 'time')
RIV_cond_winter = RIV_cond.sel(time = slice(None, "20120301")).mean(dim = 'time')

depth_spring = depths.sel(time = slice("20120301", None)).mean(dim = 'time')
RIV_stage_spring = RIV_stage.sel(time = slice("20120301", None)).mean(dim = 'time')
drainage_spring = drainage.sel(time = slice("20120301", None)).mean(dim = 'time')
riv_density_spring = riv_density.sel(time = slice("20120301", None)).mean(dim = 'time')
sy_baseline_spring = sy_baseline.sel(time = slice("20120301", None)).mean(dim = 'time')
RIV_cond_spring = RIV_cond.sel(time = slice("20120301", None)).mean(dim = 'time')

# Use gaussian blur to represent the effect of neighboring features
ds_winter = xr.Dataset({
	"transmissivity":KD,
	"resistance":VCW,
	"depth":depth_winter,
	"sto_baseline":sy_baseline_winter,
	"sto_unsat":sy,
	"sto_por":por,
	"river_stage":RIV_stage_winter,
	"river_stage_ghg":RIV_stage.isel(time = ind_max, drop = True).where(riv_flux),
	"river_stage_zeros":RIV_stage_winter.fillna(0),
	"river_stage_ghg_zeros":RIV_stage.isel(time = ind_max,drop = True).where(riv_flux).fillna(0),
	"river_stage_depth":RIV_stage_winter.fillna(depth_winter),
	"river_stage_maxdepth":RIV_stage_winter.fillna(depth_winter.max()),
	"river_conductance_ghg":RIV_cond.isel(time = ind_max, drop = True).where(riv_flux),
	"river_conductance":RIV_cond_winter,
	"river_conductance_zeros":RIV_cond_winter.fillna(0),
	"drain_binary":drainage_winter.notnull()*1.0,
	"drain_level":drainage_winter,
	"drain_level_maxdepth":drainage_winter.fillna(depth_winter.max()),
	"drain_level_depth":drainage_winter.fillna(depth_winter),
	"drain_conductivity":drainage_cond,
	"drain_conductivity_zeros":drainage_cond.fillna(0),
	"river_density":riv_density_winter,
	"river_density_ghg":riv_flux
	})
ds_winter = ds_winter.load()

ds_spring = xr.Dataset({
	"transmissivity":KD,
	"resistance":VCW,
	"depth":depth_spring,
	"sto_baseline":sy_baseline_spring,
	"sto_unsat":sy,
	"sto_por":por,
	"river_stage":RIV_stage_spring,
	"river_stage_ghg":RIV_stage.isel(time = ind_max, drop = True).where(riv_flux),
	"river_stage_zeros":RIV_stage_spring.fillna(0),
	"river_stage_ghg_zeros":RIV_stage.isel(time = ind_max,drop = True).where(riv_flux).fillna(0),
	"river_stage_depth":RIV_stage_spring.fillna(depth_spring),
	"river_stage_maxdepth":RIV_stage_spring.fillna(depth_spring.max()),
	"river_conductance_ghg":RIV_cond.isel(time = ind_max, drop = True).where(riv_flux),
	"river_conductance":RIV_cond_spring,
	"river_conductance_zeros":RIV_cond_spring.fillna(0),
	"drain_binary":drainage_spring.notnull()*1.0,
	"drain_level":drainage_spring,
	"drain_level_maxdepth":drainage_spring.fillna(depth_spring.max()),
	"drain_level_depth":drainage_spring.fillna(depth_spring),
	"drain_conductivity":drainage_cond,
	"drain_conductivity_zeros":drainage_cond.fillna(0),
	"river_density":riv_density_spring,
	"river_density_ghg":riv_flux
	})
ds_spring = ds_spring.load()

# imod.idf.save(model_dat/"min_drn_active.idf", drain_valid)
def gaussian_2d_uf(np_arr, sigma, **kwargs):
	
	if not np.isnan(np_arr).any():
		return gaussian_filter(np_arr, sigma=sigma, **kwargs)
	
	# Values of the gaussian smoothed data
	# 0 for nan and value for value
	V = np_arr.copy()
	V[np.isnan(np_arr)] = 0
	V = gaussian_filter(V, sigma=sigma, **kwargs)
	# Sum of Weights of the cells with data
	# 1s with data and 0 for nan
	W = np.ones_like(np_arr)
	W[np.isnan(np_arr)] = 0
	W = gaussian_filter(W, sigma=sigma, **kwargs)
	W[W==0]=np.nan
	
	return V/W

def gaussian_2d(da, sigma):
	n_dims = len(da.shape)
	if n_dims == 2:
		return xr.apply_ufunc(gaussian_2d_uf,
						      da, 
							  kwargs=dict(sigma = sigma),
# 							  dask="parallelized",
							  output_dtypes=[float])
	return xr.apply_ufunc(gaussian_2d_uf,
					      da,
						  kwargs = dict(sigma = sigma,
									    axes = (1,2)),
# 						  dask="parallelized",
						  output_dtypes=[float])

#%%
@delayed
def write_ds(sigma, data_set = ds_winter, csv_name = "site_properties_winter_gb_"):
	print(f"working on {sigma}")

	ds_blur = data_set.map(gaussian_2d, sigma=sigma/25)
	
	ds_blur["gaussian_sigma"] = sigma
	ds_blur = ds_blur.set_coords("gaussian_sigma")
	csv = ds_blur.to_dask_dataframe()
	
	res_dir = Path(r"C:/Users/valdf/MODFLOW_transient_res")
	csv.to_csv(res_dir/(csv_name+str(sigma)+".csv"),
				single_file=True,
				compute=True
				)
# [write_ds(n, data_set=ds).compute() for n in [1500,1000,500,250]]

delayed_compute = []
for s in [1000,500,250]:
 	delayed_compute.append(write_ds(s, data_set = ds_winter, csv_name = "site_properties_winter_gb_"))

for s in [1000,500,250]:
 	delayed_compute.append(write_ds(s, data_set = ds_spring, csv_name = "site_properties_spring_gb_"))

compute(delayed_compute, num_workers=8)
