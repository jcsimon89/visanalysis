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
# - pass roi_mask, roi_image to fly_final.hdf5 in select_final_rois - DONE
# - analyze_data final - DONE
# - save single frame pdfs and pngs
# - show stimulus times (shaded)
# - add plots: all roi plot for flash series - DONE
# - add plots: all_roi/individual roi single page summary pdf (from single panels)
# - add plots: roi image
# - SHARED_ANALYSIS.PY!!!
# - future: powerpoint slides, standardize plotting tools in analyze_data, add fano factor to gui

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
experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS140_x_JS252/fly_021' #string to folder containing fly.hdf5 file
rig = 'Bruker' #string "Bruker" or "AODscope"

roi_set_name = 'roi_set_name' # name of roi group to be analyzed (default 'roi_set_name')
response_set_name = 'mask' # name of response group to be analyzed (default 'mask')

# process_data
series_number_for_roi_selection = '1' #string 
run_gui = 'True' #string "True" or "False", default = "False"
fano_regularizer = '300' #string, k in the GUI Fano image var / (mean + k), in image intensity units
attach_metadata = 'True' #string "True" or "False", default = "False"

# analyze_data
show_figs = 'False' #string "True" or "False", default = "False"
save_figs = 'True' #string "True" or "False", default = "False"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

# select_final_rois
save_hdf5 = 'True'#string "True" or "False", default = "False"


#%% PROCESS_DATA
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell


process_data_path = str(os.path.join(base_path,'process_data.py'))

run_cmd('python ' + process_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --series_number ' + series_number_for_roi_selection
                + ' --run_gui ' + run_gui
                + ' --fano_regularizer ' + fano_regularizer
                + ' --attach_metadata ' + attach_metadata
                + ' --roi_set_name ' + roi_set_name
                + ' --response_set_name_prefix ' + response_set_name)


#%% ANALYZE_DATA RAW
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell


tag = 'raw' #string "raw" or "final"

analyze_data_path = str(os.path.join(base_path,'analyze_data.py'))

run_cmd('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff)


#%% SELECT_FINAL_ROIS
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell
input_tag = ''
output_tag = 'final'
select_rois_path = str(os.path.join(base_path,'select_rois.py'))

run_cmd('python ' + select_rois_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --save ' + save_hdf5
                + ' --input_tag ' + input_tag
                + ' --output_tag ' + output_tag)


#%% ANALYZE_DATA FINAL
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell


tag = 'final' #string "raw" or "final"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

analyze_data_path = str(os.path.join(base_path,'analyze_data.py'))

run_cmd('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff)
# %%

