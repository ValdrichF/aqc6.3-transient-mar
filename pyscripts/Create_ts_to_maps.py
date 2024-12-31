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
from matplotlib import pyplot as plt
import contextily as cx
from joblib import Parallel, delayed
from TMP_calc_response import save_response


from pathlib import Path

model_path = r"D:\Modeldatabase\Baaksebeek\MODELGRENZEN"
model_path = Path(model_path)
ibound = imod.idf.open(model_path / "MODELGRENZEN.IDF")

x = ibound.coords["x"]
y = ibound.coords["y"]

# import the original recharge file to add to the newly created one
# org_rch = imod.idf.open(model_path / "GRONDWATERAANVULLING//gemiddelde_grondwateraanvulling_new_25.IDF")

sample = pd.read_csv(r"K:/MODFLOW_transient_res/response_decay.csv")
sample["time"] = pd.to_datetime(sample.time)
times = sample.time.unique()

save_dir = (model_path.parent / "Res_maps")
save_max = save_dir/ "decay_maximumn"
save_total = save_dir/ "decay_total"
save_area = save_dir/ "decay_area"
(save_max / "pngs").mkdir(exist_ok=True, parents=True)
(save_total / "pngs").mkdir(exist_ok=True, parents=True)
(save_area / "pngs").mkdir(exist_ok=True, parents=True)

# Pretty maps
source = cx.providers.CartoDB.Voyager
kwargs_basemap = {"alpha":0.7, "zoom":10, "attribution_size":7, "attribution":""}

# Legend
legend = imod.visualize.spatial.read_imod_legend(r"D:\Modeldatabase\diff_legend.leg")
col2 = ['#03045E', '#023E8A', '#0077B6','#0096C7','#00B4D8', '#48CAE4','#90E0EF', '#ADE8F4', '#CAF0F8']
label_max = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
label_total = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
label_area = [0, 0.1, 0.25, 0.5, 0.75, 1, 2, 4]


# Loop thorugh each group
def plot_stress_period(t):
	sample_grp = sample.loc[sample['time']==t] # _no
	# Create a dummy rch file with 0s that will be filled in later
    ## Created like one of the idf files
	total = xr.full_like(ibound, np.nan, dtype = np.double)
	maximum = xr.full_like(ibound, np.nan, dtype = np.double)
	area = xr.full_like(ibound, np.nan, dtype = np.double)
	# Then through each rch_no
	for row in sample_grp.itertuples(): #_rch
		half_side = row.side*row.recharge/25/2
		selection =(
			(y   >= row.y - half_side)
			& (y   <= row.y + half_side)
			& (x   >= row.x - half_side)
			& (x   <= row.x + half_side) # sample_rch.iloc[0]
			)
		total = xr.where(selection, row.total, total) # otherwise use sample_rch.iloc[0].recharge
		maximum = xr.where(selection, row.maximum, maximum) # otherwise use sample_rch.iloc[0].recharge
		area = xr.where(selection, row.area, area) # otherwise use sample_rch.iloc[0].recharge
		
	imod.idf.save(save_total / f"total_{t.strftime('%Y%m%d')}.idf", total) # _{rech} # +org_rch
	imod.idf.save(save_max / f"maximum_{t.strftime('%Y%m%d')}.idf", maximum) # _{rech} # +org_rch
	imod.idf.save(save_area / f"area_{t.strftime('%Y%m%d')}.idf", area) # _{rech} # +org_rch
	
	fig, ax = imod.visualize.spatial.plot_map(maximum, col2[::-1], label_max, basemap=source, kwargs_basemap=kwargs_basemap)
# 	plt.title(f"Maximum $[m]$ on {t.strftime('%b-%Y')}", loc = 'right')
	ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
	plt.title("Decay of Maximum response", loc = 'right')
	plt.savefig(save_max / "pngs" / f"maximum_{t.strftime('%Y%m%d')}.png", dpi = 300)
	
	fig, ax = imod.visualize.spatial.plot_map(area, col2[::-1], label_area, basemap=source, kwargs_basemap=kwargs_basemap)
# 	plt.title(r"Area $[m^2]$ on "+t.strftime('%b-%Y'), loc = 'right')
	ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
	plt.title("Decay of Area of the response", loc = 'right')
	plt.savefig(save_area / "pngs" / f"area_{t.strftime('%Y%m%d')}.png", dpi = 300)
	
	fig, ax = imod.visualize.spatial.plot_map(total, col2[::-1], label_total, basemap=source, kwargs_basemap=kwargs_basemap)
# 	plt.title(r"Total $[m^3]$ on "+t.strftime('%b-%Y'), loc = 'right')
	ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
	plt.title("Decay of Total response", loc = 'right')
	plt.savefig(save_total / "pngs" / f"total_{t.strftime('%Y%m%d')}.png", dpi = 300)

        
plot_stress_period(times[0])

Parallel(n_jobs=12)(delayed(plot_stress_period)(n) for n in times)

plotter = save_response(Path(r"K:\MODFLOW_transient_res\avg_baseline"))

plotter.png_to_gif(save_max / "pngs")
plotter.png_to_gif(save_area / "pngs")
plotter.png_to_gif(save_total / "pngs")
        
        
        
