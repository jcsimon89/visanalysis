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
#individual rois
#experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS252/fly_012' #string to folder containing fly.hdf5 file
#with whole slice roi
#experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS252/fly_013'
experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS252/fly_014'
#experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS257/fly_007'
#experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS257/fly_008'
#experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS257/fly_009'
#experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS258/fly_001'
#experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS258/fly_002'
#experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS258/fly_003'
rig = 'Bruker' #string "Bruker" or "AODscope"

# process_data
run_gui = 'True' #string "True" or "False", default = "False"
attach_metadata = 'True' #string "True" or "False", default = "False"
show_figs = 'False' #string "True" or "False", default = "False"
save_figs = 'True' #string "True" or "False", default = "False"

# series numbers
#ex: '0,1 2,3' #series numbers for each experiment set separated by commas with no spaces, sets of series separated by space
series_num = '0,1 2,3'
#channel_num = '1,2 1,2' #channel numbers corresponding to series numbers in same format
channel_num = '1,2 1,2'
# note: first series and chennel listed will be used for roi selection

#%% PROCESS_DATA


process_data_path = str(os.path.join(base_path,'process_data_f0.py'))

os.system('python ' + process_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --run_gui ' + run_gui
                + ' --attach_metadata ' + attach_metadata
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --series_num ' + series_num
                + ' --channel_num ' + channel_num)


