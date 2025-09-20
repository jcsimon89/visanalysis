"""
visanalysis shell script for running:
1. process_data.py
2. analyze_data.py (raw)
3. select_final_rois.py
4. analyze_data.py (final)



Experiment description:
series 1: search stimulus do identify rois that are looking at the right part of the screen
series 2: rf mapping series where center location and circle radius are varied

Analysis workflow:
1. process data
2. analyze search stim data and identify good rois
3. select correct center location for each roi (if there is one, toss rois where this is not true)
4. analyze data for each roi at correct center location

using
https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu, modified by jcsimon@stanford.edu
"""
#TODO: 
# - add plots: all_roi/individual roi single page summary pdf (from single panels)
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
experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS144_x_JS252/fly_002' #string to folder containing fly.hdf5 file
rig = 'Bruker' #string "Bruker" or "AODscope"

# process_data
series_number_for_roi_selection = '1' #string 
run_gui = 'True' #string "True" or "False", default = "False"
attach_metadata = 'False' #string "True" or "False", default = "False"
roi_set_name = 'LobulaPlate' #Lobula or LobulaPlate for T5
#roi_set_name = 'LobulaPlate'
response_set_name_prefix = roi_set_name

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
                + ' --attach_metadata ' + attach_metadata
                + ' --roi_set_name ' + roi_set_name
                + ' --response_set_name_prefix '+ response_set_name_prefix)


#%% ANALYZE_DATA RAW


tag = 'raw' #string "raw" or "final"
analyze_data_path = str(os.path.join(base_path,'analyze_data_DS.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --response_set_name_prefix '+ response_set_name_prefix)

#%% SELECT_FINAL_ROIS
input_tag = ''
output_tag = 'final'
select_rois_path = str(os.path.join(base_path,'select_rois.py'))

os.system('python ' + select_rois_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --save ' + save_hdf5
                + ' --input_tag ' + input_tag
                + ' --output_tag ' + output_tag
                + ' --response_set_name_prefix '+ response_set_name_prefix)


#%% ANALYZE_DATA Final


tag = 'final' #string "raw","good","final"
analyze_data_path = str(os.path.join(base_path,'analyze_data_DS.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --response_set_name_prefix '+ response_set_name_prefix)

