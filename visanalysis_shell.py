"""
visanalysis shell script for running:
1. process_data.py
2. analyze_data.py (raw)
3. select_final_rois.py
4. analyze_data.py (final)

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
#TODO:
# 1. test analyze_data final
# 2. save single frame figs? add plots: X all roi plot for flash series, all_roi/individual roi single page summary pdf, 
# 3. future: save pngs and powerpoint slides? add roi image plots, standardize plotting tools in analyze_data, add fano factor to gui
# COMBINED ANALYSIS!!!

#%% Initialize environment

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

#%% initialize arguments

# all scripts
base_path = 'C:/Users/jcsimon/Documents/GitHub/visanalysis'
experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS139_x_JS251/fly_004' #string to folder containing fly.hdf5 file
rig = 'Bruker' #string "Bruker" or "AODscope"

# process_data
series_number_for_roi_selection = '3' #string 
run_gui = 'True' #string "True" or "False", default = "False"
attach_metadata = 'True' #string "True" or "False", default = "False"

# analyze_data
show_figs = 'False' #string "True" or "False", default = "False"
save_figs = 'True' #string "True" or "False", default = "False"

# select_final_rois
save_hdf5 = 'True'#string "True" or "False", default = "False"


#%% process_data


process_data_path = str(os.path.join(base_path,'process_data.py'))

os.system('python ' + process_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --series_number ' + series_number_for_roi_selection
                + ' --run_gui ' + run_gui
                + ' --attach_metadata ' + attach_metadata)


#%% analyze_data raw


tag = 'raw' #string "raw" or "final"
analyze_data_path = str(os.path.join(base_path,'analyze_data.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag)


#%% select_final_rois


select_final_rois_path = str(os.path.join(base_path,'select_final_rois.py'))

os.system('python ' + select_final_rois_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --save ' + save_hdf5)


#%% analyze_data final


tag = 'final' #string "raw" or "final"
analyze_data_path = str(os.path.join(base_path,'analyze_data.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag)