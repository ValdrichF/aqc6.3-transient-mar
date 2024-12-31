# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:16:22 2023

@author: Valdrich

output - csv of fluxes for each scenario
All the fluxes are in <scenario_name>/(?<metaswap>)/<FLUX_NAME>/<FLUX_NAME>_<date>_L1.IDF files
the columns would need to be:
	scenario_name
	flux_name
	date
	value

"""
from pathlib import Path
import imod
import pandas as pd
from joblib import Parallel, delayed


res_dir = Path(r"E:\MODFLOW_transient_res_wel")
big_grid = pd.read_csv(r"D:\Dropbox\WUR\Rscripts\paper_2\ortho_grid_ran_azure.csv") 

def create_ts_flux(flux_dir, scenario):
	print(f"Flux: {flux_dir.name}")
	
	baseline_dir = Path(r"E:\MODFLOW_transient_res\avg_baseline")
	
	idf = imod.idf.open(flux_dir/"*L1.idf")
	baseline = imod.idf.open(baseline_dir/"/".join(flux_dir.parts[3:])/"*idf")
	idf = idf-baseline
	idf = idf.sel(layer=1, drop=True)
	
	# Only where there's a response >1cm
	scenario_dir = list(flux_dir.parents)[-3]
	
	nearest_rch = imod.idf.open(scenario_dir/"nearest_rch"/"rch_l1.idf")
	
	scenario_num = int(scenario_dir.name.split("_")[1])
	
	LHS_subset = big_grid.loc[(big_grid["scenario"] == scenario_num)]
	
	table = []
	for i, row in enumerate(LHS_subset.itertuples()):
		idf_interest = idf.where(nearest_rch==row.recharge)
		total_flux = idf_interest.sum(dim=["x", "y"])
		
		total_flux_df = total_flux.to_dataframe(name="value")
		
		total_flux_df["flux_name"] = flux_dir.name
		total_flux_df["response_area"] = (nearest_rch==row.recharge).sum().values*25*25/1000/1000
		total_flux_df["scenario"] = scenario
		total_flux_df["rch_no"] = i
		total_flux_df["recharge"] = row.recharge
		total_flux_df["area"] = row.area
		table.append(total_flux_df)
	
	return pd.concat(table)


def create_ts_scenario(scenario_dir):
	
	if False:#(scenario_dir/"total_fluxes.csv").exists():
		print(f"{scenario_dir/'total_fluxes.csv'} already exists. Skipping...")
		return
	
	print(f"Working on scenario: {scenario_dir.name}")
	
	results_dirs = [n.parent for n in scenario_dir.rglob("*L1.idf")]
	results_dirs = [n for n in results_dirs if (n.name!="nearest_rch")&(n.name!="response")]
	results_dirs = list(set(results_dirs))

	all_df = [create_ts_flux(df, scenario=scenario_dir.name) for df in results_dirs]
	
	all_df = pd.concat(all_df)
	
	all_df.to_csv(scenario_dir/"total_fluxes_sites.csv")

all_scenarios = [n for n in res_dir.iterdir() if n.is_dir() & (n.name!="avg_baseline")]

Parallel(n_jobs=5)(delayed(create_ts_scenario)(n) for n in all_scenarios)








