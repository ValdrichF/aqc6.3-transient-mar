# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:33:43 2023

@author: Valdrich

output - csv of avg KD, VCW, and the depth to groundwater within 1km from the recharge site

	Drain and Rivers, only consider those that are active, 
	ie stage below the groundwater level or infiltration factor == 1
	
"""
from pathlib import Path
import imod
import pandas as pd
import xarray as xr
from nearest_sites import nearest_site

res_dir = Path(r"Z:\MODFLOW_transient_res")
big_grid = pd.read_csv(r"D:\Dropbox\WUR\Rscripts\paper_2\ortho_grid.csv") 

# maximum baseline head
baseline_heads = imod.idf.open(res_dir/"avg_baseline/head/*l1.idf")
baseline_heads = baseline_heads.sel(layer = 1, drop=True)

# KD, VCW
small_model = Path(r"D:\Modeldatabase\Baaksebeek\Small_model")
KD = imod.idf.open(small_model/"Transmissivity/Transmissivity_l1.idf")
VCW = imod.idf.open(small_model/"Resistance/Resistance_l1.idf")
KD = KD.sel(layer = 1, drop = True)
VCW = VCW.sel(layer = 1, drop = True)

# Depths
depths = imod.idf.open(res_dir/"avg_baseline/depth/*l1.idf")
depths = depths.sel(layer = 1, drop=True)

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

Xmin = baseline_heads.x.min().values
Xmax = baseline_heads.x.max().values
Ymin = baseline_heads.y.min().values
Ymax = baseline_heads.y.max().values

RIV_stage = RIV_stage.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin))
RIV_cond = RIV_cond.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin))
drainage = drainage.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin))
drainage_cond = drainage_cond.sel(x = slice(Xmin,Xmax), y = slice(Ymax, Ymin))

# approximate river density, normally drain density, but AMIGO
riv_density = RIV_stage.notnull()

# Select the cells with a river in that stress-period
RIV_cond = RIV_cond.where(riv_density)

# Filling nan conductace with 0s
RIV_cond = RIV_cond.fillna(0)
drainage_cond = drainage_cond.fillna(0)

# imod.idf.save(model_dat/"stage_active.idf", RIV_stage_valid)

# imod.idf.save(model_dat/"min_drn_active.idf", drain_valid)
ds = KD.to_dataset(name = "transmissivity")
ds["resistance"] = VCW
ds["depth"] = depths
ds["river_stage"] = RIV_stage
ds["river_conductance"] = RIV_cond
ds["drain_level"] = drainage
ds["drain_conductivity"] = drainage_cond
ds["river_density"] = riv_density
ds.load()

# Distance from the recharge site
n_site = nearest_site()

table = []
# Extract the avg within 1km from the site
for s in big_grid.scenario.unique():
	LHS_subset = big_grid.loc[(big_grid["scenario"] == s)]

	for row in LHS_subset.itertuples():
		print(f"working on scenario: {row.scenario}, rch_site: {row.rch_site}")
		dist_site = n_site.dist_to_rch(row.N_side, row.S_side, row.W_side, row.E_side)
		
		# Static properties (only dependent on distance & time)
		arr = ds.where(dist_site<=1000).mean(dim = ["x", "y"])
		
		csv = arr.to_dataframe()
		csv["scenario"] = row.scenario
		csv["rch_site"] = row.rch_site
		table.append(csv)

geo_properties = pd.concat(table)
geo_properties.to_csv(res_dir/"site_properties_RivCondPerPeriod.csv")
