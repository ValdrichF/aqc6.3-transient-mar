# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 16:36:23 2024

@author: Valdrich

Plot a column from a csv like big_grid

Need to plot 
	- the most efficient recharge rate
	- the AMIGO recharge rate
	- the maximum ML efficiency
	- the efficiency of the AMIGO scenario

"""

# import the packages
import imod 
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colormaps as cmaps
import contextily as cx
from joblib import Parallel, delayed


from pathlib import Path

model_path = r"E:\Modeldatabase\Baaksebeek\MODELGRENZEN"
model_path = Path(model_path)
ibound = xr.open_dataset(r"E:\Modeldatabase\data_scenarios_New2\train1000\x\rch_1_1_0.nc")
ibound = xr.ones_like(ibound.rch)
x = ibound.coords["x"]
y = ibound.coords["y"]

# import the original recharge file to add to the newly created one
# org_rch = imod.idf.open(model_path / "GRONDWATERAANVULLING//gemiddelde_grondwateraanvulling_new_25.IDF")

sample_dir = Path(r"Z:/stationary_diff_rch_10ha")
csv_name = "rch_80pc_percent"
sample = pd.read_csv(sample_dir/f"{csv_name}.csv")

(sample_dir/"idf").mkdir(exist_ok = True)

# Pretty maps
source = cx.providers.CartoDB.Voyager
kwargs_basemap = {"alpha":1, "zoom":10, "attribution_size":7, "attribution":""}

# Legend
legend = imod.visualize.spatial.read_imod_legend(r"E:\Modeldatabase\Baaksebeek\diff_legend.leg")
col2 = ['#03045E', '#023E8A', '#0077B6','#0096C7','#00B4D8', '#48CAE4','#90E0EF', '#ADE8F4', '#CAF0F8']
label = [60, 120, 180, 365, 1.5*365.25, 2*365, 3*365+1, 5*365+1]
col_rch = ['#f1eef6', '#d0d1e6', '#a6bddb', '#74a9cf', '#2b8cbe', '#045a8d']
label_rch = list(range(5, 26, 2))

cmap = cmaps['winter']
col_rch = cmap(np.linspace(0, 1, len(label_rch)+1))
col_rch = col_rch[::-1]



# Loop thorugh each group
def plot_col(column_name, colours, legend, dat = sample):
	# Create a dummy rch file with 0s that will be filled in later
    ## Created like one of the idf files
	total = xr.full_like(ibound, np.nan, dtype = np.double)
	# Then through each rch_no
	for row in dat.itertuples(): #_rch
		half_side = row.area_site**0.5/2
		selection =(
			(y   >= row.y - half_side)
			& (y   <= row.y + half_side)
			& (x   >= row.x - half_side)
			& (x   <= row.x + half_side) # sample_rch.iloc[0]
			)
		total = xr.where(selection, getattr(row, column_name), total) # otherwise use sample_rch.iloc[0].recharge
		
	imod.idf.save(sample_dir/"idf" / f"{csv_name}_{column_name}.idf", total) # _{rech} # +org_rch
		
	fig, ax = imod.visualize.spatial.plot_map(total, colours, legend, basemap=source, kwargs_basemap=kwargs_basemap)
# 	plt.title(r"Total $[m^3]$ on "+t.strftime('%b-%Y'), loc = 'right')
	ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
	plt.title("", loc = 'right')
	plt.savefig(sample_dir / f"{csv_name}_{column_name}.png", dpi = 300)
	plt.clf()
	


plot_col('rch', col_rch, label_rch, sample)
plot_col('ideal_rch', col_rch, label_rch, sample)


col2 = ['#CAF0F8', '#023E8A', '#0077B6','#00B4D8', '#48CAE4', '#ADE8F4']
label2 = [5, 10, 15, 20, 25]
fig, ax = imod.visualize.spatial.plot_map(total, col2[::-1], label2, basemap=source, kwargs_basemap=kwargs_basemap)
ax.tick_params(which='both', labelbottom=False, labelleft=False) #  bottom=False, left=False, 
plt.title("", loc = 'right')
plt.savefig(sample_dir / "rch_80pc_percent_rch2.png", dpi = 1000)
plt.clf()
