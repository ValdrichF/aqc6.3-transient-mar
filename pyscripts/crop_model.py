# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 11:55:49 2022

@author: Valdrich

Crop all the input files to the Baaksebeek catchment area +/- spreading length
The spreading length for the phreatic aquifer is collected from: Altera report 1102
https://www.researchgate.net/publication/40796941_Monitoring_van_verdroging_methodische_aspecten_van_meetnetoptimalisatie

for now, only the 25m stationary data are cropped
"""

import imod
from pathlib import Path
from os import makedirs
import xarray as xr
from shutil import copyfile, copytree
import re

model_path = r'D:\Modeldatabase\AMIGO31'
model_path = Path(model_path)

cropped_path = r'E:\Modeldatabase\AMIGO_needed'
cropped_path = Path(cropped_path)

if not cropped_path.is_dir():
    makedirs(cropped_path)

sample = imod.idf.open(r"D:\Modeldatabase\Baaksebeek\MODELGRENZEN\MODELGRENZEN.IDF")
"""
Transmissivity, Vertical Resistance, SHD, Horizontal Anisotropy, Recharge, 
Drainage, Overland flow, Rivers

Wells is an .ipf file

Loop through the folders, loop through each idf in it, crop it and then save it
"""
# Loop through each folder
# folders = ["ANISOTROPIE", "C-WAARDEN", "DRAINAGE", "KD-WAARDEN", "MEETREEKSEN",
# 		   "METASWAP", "MODELGRENZEN", "OLF", "ONTTREKKINGEN", "OPPERVLAKTEWATER",
# 		   "rch_transient", "MAAIVELD", "CHD", "METEO", "METASWAP\Initialisatie"]
folders = ["OPPERVLAKTEWATER"]
## find the path to all the idf files within the folder
idfs = model_path.rglob("*")

## Loop through each idf 
for folder in folders:
	print(folder)
	
	for file in (model_path/folder).rglob("*"):
		
		if file.is_dir():
			continue
		new_folder = cropped_path / "/".join(file.parts[3:-1])
		new_folder.mkdir(parents=True, exist_ok=True)
		save_path = new_folder / file.parts[-1]
		
		print(file)

		suffix = file.suffix.lower()
		if (suffix == ".zip")|(suffix == ".asc"):
			continue
		if (suffix==".idf"):
			
			### Read the idf
			data = imod.idf.open(file)
			### Remove the unnecessary dimensions so that imod doesn't change the file name
			data = data.squeeze()
			if "layer" in data.coords._names:
				data = data.drop_vars(["layer"])
			if "time" in data.coords._names:
				data = data.drop_vars(["time"])
			
			# Reindex and interpolate like the sample
			data_cropped = data.interp_like(sample, method='nearest')
			
			### Save it
			imod.idf.write(save_path, data_cropped)
			
			continue
		
		copyfile(file, save_path)

# Copy CHD and meteo folders as they have a different resolution.
# interpolating the values increases the size.
# Relying on iMOD's model cropping
folders = ["CHD", "METEO"]
for folder in folders:
	print(folder)
	
	copytree(model_path/folder, cropped_path/folder)

#  (y: 1009, x: 1647)
"""
 Boundary
 
 load one of the data files (recharge)
 create a new idf with only ones with the same dims as prepared file
 edit the edges to -1 ie. constant head 
"""

bnd = imod.idf.open(cropped_path/"MODELGRENZEN"/"MODELGRENZEN.IDF")

bnd = xr.ones_like(bnd)

bnd[0] = -1
bnd[-1] = -1
bnd[:,0] = -1
bnd[:,-1] = -1

imod.idf.save(cropped_path / "MODELGRENZEN" / "MODELGRENZEN.IDF", bnd)