# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:43:12 2022

@author: Valdrich

Use the created model to predict the change in the head of the test set.
These results will be saved in a folder with a similar name to the model_dir name
"""
import imod
import keras
from pathlib import Path
from functions import load_img_separate
import xarray as xr
import matplotlib.pyplot as plt
import contextily as cx
import numpy as np
import pandas as pd
import plotnine as p9
from timeit import default_timer as timer

#%% Define some models to help plot the results

def plot_results(df, model_path, title = "", x_lab = "Total recharge (m3/day)", y_lab = ""):

	df["total_recharge"] = df.mean_rch * df.area * 1000
	df = pd.wide_to_long(df,stubnames="stub", i = ["test_name", "total_recharge"], j="type", sep = "_", suffix="\\w+")
	df = df.reset_index(level=[0,1,2])
	plot = (p9.ggplot(df)+
		p9.geom_point(p9.aes(x = "total_recharge", y = "stub", color = "type"))+
		p9.labs(x = x_lab, y = y_lab, title = title)+
		p9.theme_light()+
		p9.theme(legend_title=p9.element_blank())
		)
	plot.save(model_path / (title + ".png"), dpi = 300)

def plot_pred(df, model_path, title = "AMIGO vs ML model", x_lab = "AMIGO", y_lab = "Proposed model"):
	plot = (p9.ggplot(df)+
		 p9.geom_point(p9.aes(x = "stub_Actual", y = "stub_Predicted"))+
		 p9.labs(x = x_lab, y = y_lab, title = title)+
		 p9.geom_abline(p9.aes(slope = 1, intercept = 0,linetype="'a'"))+
		 p9.theme_light()+
		 p9.theme(legend_title=p9.element_blank()) +
		 p9.scale_linetype_manual(values = dict(a = "dotted"), labels = "one to one")
		 )
	plot.save(model_path / (title + "VS.png"), dpi = 300)


#%% Loop through each model_folders in model_dir
model_dir = Path(r"D:\Modeldatabase\Azure_output_Compiled")

# List out all the directories (each model in its own dir)
model_folders = [x.parent for x in model_dir.rglob("history.csv")]

# test_dir = Path(r"D:\Modeldatabase\Haaksebergen\rch_test")
# test_img_list = [n for n in (test_dir/"x").glob("*.nc")] #  / "x"

# diff = imod.idf.open(r"D:\Modeldatabase\Haaksebergen\rch_test\y\rch_1.idf")

# #
# nc_path = Path(r"D:\Modeldatabase\Haaksebergen\common_ncs_New")
# common = xr.open_dataset(nc_path / "common_data_DE_KD_VR_DC_DL.nc").to_array()
# common = common.transpose(..., "variable")
# common = common.values.reshape((1, *common.shape))

# test_dir = Path(r"D:\Modeldatabase\data_scenarios_New2\test")
# test_dir = Path(r"F:\Modeldatabase\AMIGO_needed\rch_transient_individual_sites\nc")
test_dir = Path(r"F:\Modeldatabase\AMIGO_needed\rch_transient\cropped")
test_img_list = list(test_dir.glob("*.nc")) #  / "x"

# diff = xr.open_dataset(test_dir/"y/rch_5_3_4.nc")

#
nc_path = Path(r"D:\Modeldatabase\common_ncs_New")
common = xr.open_dataset(nc_path / "common_data_DE_KD_VR_DC_DL.nc").to_array()
common = common.transpose(..., "variable")
# common = common.values.reshape((1, *common.shape))



# len_x = len(diff.x)
# len_y = len(diff.y)

# ## To scale the depth to groundwater back from 0-1 to min_diff to max_diff
# min_diff = 0
# max_diff = 10
# log_diff = False

# Some constants for plotting
Xmin = 204519
Xmax = 241200

Ymin = 439298
Ymax = 460020

# Background
source = cx.providers.CartoDB.Voyager
kwargs_basemap = {"alpha":0.7, "zoom":10, "attribution_size":2, "attribution":""}

# Legend
legend = imod.visualize.spatial.read_imod_legend(r"D:\Modeldatabase\diff_legend.leg")
error_col = ['#8c510a','#d8b365','#f6e8c3','#ffffff','#c7eae5','#5ab4ac','#01665e']
error_label = [-0.5, -0.1, -0.01, 0.01, 0.1, 0.5]
ibound = common.sel(variable = "depth", drop = True)
ibound = xr.ones_like(ibound)

model_folder = model_folders[5]
    # use  model_folder.parent.name[6:]) if the model is trained locally
    # Use list(model_folder.parents)[1].name is the model is trained using azureml
preds_dir = Path(r"Z:\Stationary_preds_multiple") / list(model_folder.parents)[1].name / f"preds_{model_folder.name}"
(preds_dir/"compare_png").mkdir(parents = True, exist_ok=True)


# Find the lastest hdf5 file in the directory
list_of_files = list(model_folder.glob('*.hdf5'))
latest_file = max(list_of_files)

model = keras.models.load_model(latest_file, compile = False)

# Some lists to save the calculated metrics
limit = 0.1 # we're interested in the extent of response greater than or equal to...
max_actual = []
max_pred = []
extent_x_actual = []
extent_y_actual = []
extent_x_pred = []
extent_y_pred = []
area_actual = []
area_pred = []
tot_actual = []
tot_pred = []
rate = []
area = []
errors = []
null_error = []

roots = []

# list to save times in
times = []

# stopper = 0
for count, i in enumerate(test_img_list):
	print(count)
    # if stopper > 4:
    #     break
    # stopper += 1

	root = i.name
	print("working on:", root)

    ## Save the predictions
	test_img = load_img_separate([i], common)
	start = timer()
	test_pred = model.predict(test_img)
	times.append(timer()-start)

	diff_save =  xr.DataArray(test_pred[0,:, :, 0], coords=ibound.coords, dims = ibound.dims)
	diff_save = diff_save.where(diff_save>0,0)

# 	imod.idf.save(preds_dir / "Predictions" / ("pred" + i.stem[3:] + ".idf"),diff_save)
	imod.idf.save(preds_dir / "Predictions" / (i.stem + ".idf"),diff_save)
    # Calculate the error to be saved as part of the png

    # Save as png only for individual recharge sites

	## Count the number of digits in the file name
	## 3 = individual sites, 4 = pairs of sites
	## n_dig == 2 because root.split doesn't separate 4.nc for eg.
	n_dig = len([num for num in root.split(sep = "_") if num.isdigit()])

	if True:#n_dig == 2:

	    fig, axes = plt.subplots(nrows=1, sharex=True, figsize = (6,4))
	    imod.visualize.spatial.plot_map(diff_save, legend[0], legend[1], basemap=source, fig = fig, ax = axes, kwargs_basemap=kwargs_basemap)
	    axes.set_title("Stationary response")
# 		    axes[2].set_title("Error\nProposed model - AMIGO")
# 		    axes[0].tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
# 		    axes[1].tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
# 		    axes[2].tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False) #  bottom=False, left=False, 
	    plt.tight_layout()
	    plt.savefig(preds_dir / "compare_png" / ("res_compare"+root[3:-3]+".png"), dpi = 300)
	    plt.draw()
	    plt.clf()

    #%% Calculate the metrics to judge the model performance
	rch  = xr.open_dataset(i)["rch"]
	rate.append(
		rch.where(rch>0).mean()
	)
	area.append(
		(rch>0).sum()*rch.dx*rch.dy*-1/1000/1000
	)
	# For the max
	max_pred.append(
		float(diff_save.max())
	)

    # Extent x
	## Fist calculate a 1D array of the maximum response along the y direction
	## Then count how many cells are higher than the threshold limit (bool sum)
	## Then convert it to km based on the grid size

	d = (diff_save.max(dim = "y")>limit).sum()
	extent_x_pred.append(
        (d*diff_save.dx/1000).values
    )

	# Extent y
	d = (diff_save.max(dim = "x")>limit).sum()
	extent_y_pred.append(
        (d*diff_save.dy/1000).values*-1
    )
	
	# Area
	d = (diff_save>0.05).sum()
	area_pred.append(
		d*(diff_save.dy/1000)**2
	)

    # total
	tot_pred.append(
		(diff_save.sum()*diff_save.dx*diff_save.dy).values*-1
	)

    # Save the name of the root (scenario) for later reference
	roots.append(root)

# Save the summary metrics from all the test data
rate = np.array(rate)
area = np.array(area)
extent_x = pd.DataFrame(dict(test_name = roots,
							 mean_rch = rate,
							 area = area,
                             stub_Predicted = np.array(extent_x_pred),
                             stub_Actual = np.array(extent_x_actual)))
extent_y = pd.DataFrame(dict(test_name = roots,
							 mean_rch = rate,
							 area = area,
                             stub_Predicted = np.array(extent_y_pred),
                             stub_Actual = np.array(extent_y_actual)))
total = pd.DataFrame(dict(test_name = roots,
						  mean_rch = rate,
						  area = area,
                          stub_Predicted = np.array(tot_pred),
                          stub_Actual = np.array(tot_actual)))
highest = pd.DataFrame(dict(test_name = roots,
							mean_rch = rate,
							area = area,
							stub_Predicted = np.array(max_pred),
							stub_Actual = np.array(max_actual)))
results = pd.DataFrame(dict(test_name = roots,
							mean_rch = rate,
							area = area,
                            stub_UNET = np.array(errors),
                            stub_Null_model = np.array(null_error)))
area_df = pd.DataFrame(dict(test_name = roots,
							mean_rch = rate,
							area = area,
                            stub_UNET = np.array(area_pred),
                            stub_Actual = np.array(area_actual)))
times = pd.DataFrame(times, columns=["time"])
times.to_csv(preds_dir / "times.csv", index = False)


plot_results(extent_x, preds_dir, "E-W extent of response", y_lab = "Length (km)")
plot_results(extent_y, preds_dir, "N-S extent of response", y_lab = "Length (km)")
plot_results(total, preds_dir, "Total response", y_lab = "Volume (m3)")
plot_results(highest, preds_dir, "Maximum response", y_lab = "Increase (m)")
plot_results(results, preds_dir, "Model error vs NULL model", y_lab = "RMSE (m)")

plot_pred(extent_x, preds_dir, "E-W extent of response")
plot_pred(extent_y, preds_dir, "N-S extent of response")
plot_pred(total, preds_dir, "Total response")
plot_pred(highest, preds_dir, "Maximum response")

#%% Special plot for error. Normalised by the maximum response
results["highest"] = highest.stub_Actual
results["norm_rmse"] = results.stub_UNET/results.highest
plot = (
        p9.ggplot(results)+
        p9.geom_point(p9.aes(x = "total_recharge", y = "norm_rmse"), color = "steelblue")+
        p9.labs(x = "Total recharge (m3/day)",
				y = "RMSE/maximum response (-)",
				title = "Normalised error")
        )
plot.save(preds_dir / "Normalised error.png", dpi=300)


extent_x.to_csv(preds_dir / "extent_x.csv", index=False)
extent_y.to_csv(preds_dir / "extent_y.csv", index=False)
total.to_csv(preds_dir / "total.csv", index=False)
highest.to_csv(preds_dir / "highest.csv", index=False)
results.to_csv(preds_dir / "rmse_error.csv", index=False)
area_df.to_csv(preds_dir / "area.csv", index=False)
