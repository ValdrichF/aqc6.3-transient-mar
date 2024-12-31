# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 14:59:28 2023

@author: Valdrich
calculate the difference between the baseline scenario and the recharge scenario
Given a directory name, do it...
"""

from pathlib import Path
import imod

scenarios_dir = Path(r"D:\MODFLOW_transient_res")

baseline = "avg_initial"
result_dir = "head"
save_dir = "differences"

baseline_idf = imod.idf.open(scenarios_dir/baseline/result_dir/"*.idf")

scenarios = [n for n in scenarios_dir.iterdir() if n.is_dir()]

def difference_scenario(path, baseline_idf, save_dir):
	
	scenario_idf = imod.idf.open(path/"*.idf")
	
	difference = scenario_idf-baseline_idf
	
	difference.rename(f"change in {path.parent.name}")
	
	save_dir.mkdir(exist_ok=True, parents=True)
	
	imod.idf.save(save_dir, difference)
	
[difference_scenario(n/result_dir, baseline_idf, scenarios_dir/save_dir/result_dir/n.name/result_dir) for n in scenarios]

