"""
visanalysis shell script for running:
1. process_data.py
2. analyze_data.py (raw)
3. select_final_rois.py
4. analyze_data.py (final)

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
#TODO: pass roi_mask, roi_image to fly_final.hdf5 in analyze_data (raw)
# plot roi_image in process_data before and after saveroiresponsefrommask (is it just a meanbrain image the whole time?? i think so)
# decide how to save mask (single array with different number for each mask or list of individual masks?)

# - test analyze_data final
# - save single frame figs
# - show stimulus times (shaded)
# - add plots: all roi plot for flash series - DONE
# - add plots: all_roi/individual roi single page summary pdf (from single panels)
# - SHARED_ANALYSIS.PY!!!
# - future: save pngs and powerpoint slides? add roi image plots, standardize plotting tools in analyze_data, add fano factor to gui

#%% INITIALIZE ENVIRONMENT

import sys
import os
import argparse
import json
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from visanalysis.plugin import base as base_plugin
from visanalysis.analysis import imaging_data
import h5py

#%% INITIALIZE ARGUMENTS

# all scripts
base_path = 'C:/Users/jcsimon/Documents/GitHub/visanalysis'
experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS139_x_JS251/fly_004' #string to folder containing fly.hdf5 file
rig = 'Bruker' #string "Bruker" or "AODscope"

# process_data
series_number_for_roi_selection = '3' #string 
run_gui = 'False' #string "True" or "False", default = "False"
attach_metadata = 'False' #string "True" or "False", default = "False"

# analyze_data
show_figs = 'False' #string "True" or "False", default = "False"
save_figs = 'True' #string "True" or "False", default = "False"

# select_final_rois
save_hdf5 = 'True'#string "True" or "False", default = "False"


#%% PROCESS_DATA


process_data_path = str(os.path.join(base_path,'process_data.py'))

os.system('python ' + process_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --series_number ' + series_number_for_roi_selection
                + ' --run_gui ' + run_gui
                + ' --attach_metadata ' + attach_metadata)


#%% ANALYZE_DATA RAW


tag = 'raw' #string "raw" or "final"
analyze_data_path = str(os.path.join(base_path,'analyze_data.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag)


#%% SELECT_FINAL_ROIS


select_final_rois_path = str(os.path.join(base_path,'select_final_rois.py'))

os.system('python ' + select_final_rois_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --save ' + save_hdf5)


#%% ANALYZE_DATA FINAL


tag = 'final' #string "raw" or "final"
analyze_data_path = str(os.path.join(base_path,'analyze_data.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag)
# %%
