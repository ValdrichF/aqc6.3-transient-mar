# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:15:08 2023

@author: Valdrich



"""

from pathlib import Path
import imod
from matplotlib import pyplot as plt
import contextily as cx
import imageio
import pandas as pd
from joblib import Parallel, delayed
import re
from pygifsicle import optimize

class save_response():

	def __init__(self, baseline_dir):
		baseline = imod.idf.open(baseline_dir/"head"/"*L1.idf")
		self.idf_baseline = baseline.sel(x = slice(204519, 241200), y = slice(460020, 439298))
		baaksebeek_model = Path(r"D:\Modeldatabase\Baaksebeek")
		self.legend = imod.visualize.spatial.read_imod_legend(baaksebeek_model / "diff_legend.leg")
# 		legend_colors = ['#ffffff00']
# 		legend_colors.extend(self.legend[0][1:])
# 		self.legend_color = legend_colors
		self.source = cx.providers.CartoDB.Voyager
		self.kwargs_basemap = {"zoom":10, "attribution_size":2, "attribution":""}
		self.duration = [.25]*21 + [1] + [.25]*30 + [1]
		

	def to_idf(self, scenario_path):
		
	
		idf_scenario = imod.idf.open(scenario_path/"head"/"*L1.idf")
		
		diff = idf_scenario-self.idf_baseline
		
# 		imod.idf.save(scenario_path/"response"/"response", diff)
		
		# save as png
		save_dir = scenario_path / "scenarios_png_L1"
		if save_dir.exists(): 
			[n.unlink() for n in save_dir.iterdir()]
		else:
			save_dir.mkdir()
		
		diff = diff.sel(layer=1, drop=True)
		
		for i in diff.time:
			i_date = i.dt.strftime('%Y%m%d').values
			i_date_title = i.dt.strftime('%m-%Y').values
			print(f"working on: {i_date}")
			diff_plot = diff.sel(time = i)
			diff_plot = diff_plot.where(diff_plot > self.legend[1][0])
			print(diff_plot.min().values)
			fig, ax = imod.visualize.spatial.plot_map(diff_plot,
													   self.legend[0], 
													   self.legend[1], 
													   basemap=self.source, 
													   kwargs_basemap=self.kwargs_basemap)
			original_x_ticks = ax.get_xticks()[1:-1]
			new_x_ticks = (original_x_ticks/1000).astype(int)
			ax.set_xticks(original_x_ticks)
			ax.set_xticklabels(new_x_ticks)

			original_y_ticks = ax.get_yticks()[1::2]
			new_y_ticks = (original_y_ticks/1000)
			ax.set_yticks(original_y_ticks)
			ax.set_yticklabels(new_y_ticks)
			
			ax.set_title(f"Response on {i_date_title}")
			plt.tight_layout()
			plt.savefig(save_dir / f"response_{i_date}.png", dpi = 300)
			plt.draw()
			plt.clf()
		
		

	def png_to_gif(self, png_dir):
		
		imgs = sorted(png_dir.iterdir(), key = lambda x:pd.to_datetime(re.findall(r'\d+', x.stem), dayfirst=False))
		
		imgs_list = []
		for f in imgs:
			image = imageio.imread(f)
			imgs_list.append(image)
		imgs_list = [imageio.imread(n) for n in imgs]
		imageio.mimsave(png_dir.parent/"ts.gif",
				        imgs_list,
						subrectangles = True,
						duration = self.duration)
		optimize(png_dir.parent/"ts.gif")

	
	def save_all(self, scenario_path):
		print(f"Starting with scenario: {scenario_path.name}")
		self.to_idf(scenario_path)
		self.png_to_gif(scenario_path/"scenarios_png_L1")
		
if __name__ == "__main__":
	baseline = Path(r"Z:\MODFLOW_transient_res\avg_baseline")
	scenarios = [n for n in baseline.parent.iterdir() if n.is_dir()]
	del scenarios[0:3]
	
	# [plotter.save_all(n) for n in scenarios]
	plotter = save_response(baseline)
	
	Parallel(n_jobs=7)(delayed(plotter.save_all)(n) for n in scenarios)





