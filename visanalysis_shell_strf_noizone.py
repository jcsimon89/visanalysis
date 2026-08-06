"""
visanalysis shell script for running:
1. process_data.py (series 1: identify rois from search/flash stim)
2. analyze_data.py (raw)
3. select_final_rois.py
4. analyze_data.py (final)
5. analyze_data_strf_noizone.py (final) - STRF from series 2 noizone stim, using final rois

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
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
experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS140_x_JS261/fly_010' #string to folder containing fly.hdf5 file
rig = 'Bruker' #string "Bruker" or "AODscope"

roi_set_name = 'roi_set_name' # name of roi group to be analyzed (default 'roi_set_name')
response_set_name = 'mask' # name of response group to be analyzed (default 'mask')

# process_data
series_number_for_roi_selection = '1' #string, first series (search/flash stim) is used to identify rois
run_gui = 'True' #string "True" or "False", default = "False"
attach_metadata = 'True' #string "True" or "False", default = "False"
flymax_movies_root = 'C:/Users/jcsimon/Documents/Stanford/Data/flymax_movies' #root dir to search for noise stimulus .bin/_params.json files

# analyze_data
show_figs = 'False' #string "True" or "False", default = "False"
save_figs = 'True' #string "True" or "False", default = "False"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

# select_final_rois
save_hdf5 = 'True'#string "True" or "False", default = "False"

# analyze_data_strf_noizone
filter_length = '5' #string, STRF filter length in seconds
monitor_hz = '' #string, monitor refresh rate override in Hz, leave '' to use ImagingDataObject default (120 Hz)


#%% PROCESS_DATA


process_data_path = str(os.path.join(base_path,'process_data.py'))

os.system('python ' + process_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --series_number ' + series_number_for_roi_selection
                + ' --run_gui ' + run_gui
                + ' --attach_metadata ' + attach_metadata
                + ' --roi_set_name ' + roi_set_name
                + ' --response_set_name_prefix ' + response_set_name
                + ' --flymax_movies_root ' + flymax_movies_root)


#%% ANALYZE_DATA RAW


tag = 'raw' #string "raw" or "final"

analyze_data_path = str(os.path.join(base_path,'analyze_data.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff)


#%% SELECT_FINAL_ROIS
input_tag = ''
output_tag = 'final'
select_rois_path = str(os.path.join(base_path,'select_rois.py'))

os.system('python ' + select_rois_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --save ' + save_hdf5
                + ' --input_tag ' + input_tag
                + ' --output_tag ' + output_tag)


#%% ANALYZE_DATA FINAL


tag = 'final' #string "raw" or "final"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

analyze_data_path = str(os.path.join(base_path,'analyze_data.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff)


#%% ANALYZE_DATA_STRF_NOIZONE FINAL
# uses the final rois selected above; auto-detects the noizone series (series 2) within the fly

tag = 'final' #string "raw" or "final"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

analyze_data_strf_path = str(os.path.join(base_path,'analyze_data_strf_noizone.py'))

command = ('python ' + analyze_data_strf_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff
                + ' --filter_length ' + filter_length)
if monitor_hz != '':
    command += ' --monitor_hz ' + monitor_hz

os.system(command)
# %%
