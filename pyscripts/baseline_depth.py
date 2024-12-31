# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 17:33:56 2023

@author: valdf

depth to the baseline groundwater head
"""

from pathlib import Path
import imod


res_dir = Path(r"K:\MODFLOW_transient_res/avg_baseline")

baseline_head = imod.idf.open(res_dir/"head/*l1.idf")
olf = imod.idf.open(Path(r"D:\Modeldatabase\Baaksebeek")/"OLF/OLF.idf")


depths = olf-baseline_head

imod.idf.save(res_dir/"depth"/"depth", depths.transpose(..., "y", "x"))
