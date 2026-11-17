"""
visanalysis shell script for running:
1. process_data.py (series 1: identify rois from search/flash stim)
2. analyze_data.py (raw)
3. select_final_rois.py
4. analyze_data.py (final)
5. analyze_data_strf_one_d.py (final) - STRF from the single-layer 1D bar/wedge
   ternary noise series, in the cell basis, using final rois

Same block structure as visanalysis_shell_strf_octave.py: run it top to bottom
the first time, then re-run individual blocks (#%% cells) while tuning.

https://github.com/ClandininLab/visanalysis
"""
#%% INITIALIZE ENVIRONMENT

import os

#%% INITIALIZE ARGUMENTS

# all scripts
base_path = 'C:/Users/jcsimon/Documents/GitHub/visanalysis'
experiment_file_directory = 'C:\\Users\\jcsimon\\Documents\\Stanford\\Data\\Bruker\\eyesss\\JS140_x_JS261\\fly_021' #string to folder containing fly.hdf5 file
rig = 'Bruker' #string "Bruker" or "AODscope"

roi_set_name = 'roi_set_name' # name of roi group to be analyzed (default 'roi_set_name')
response_set_name = 'mask' # name of response group to be analyzed (default 'mask')

# process_data
series_number_for_roi_selection = '1' #string, first series (search/flash stim) is used to identify rois
run_gui = 'True' #string "True" or "False", default = "False"
fano_regularizer = '300' #string, k in the GUI Fano image var / (mean + k), in image intensity units
attach_metadata = 'True' #string "True" or "False", default = "False"
flymax_movies_root = 'C:/Users/jcsimon/Documents/Stanford/Data/flymax_movies' #root dir to search for noise stimulus .bin/_params.json files
                                                                              #unused by the 1D stimulus, which has no movie file (it is
                                                                              #regenerated from the epoch metadata, same as the octave one),
                                                                              #but process_data still takes it for other series
monitor_hz = '' #string, display rate override in Hz, leave '' to use the ImagingDataObject default (120 Hz)


#%% PROCESS_DATA
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell


process_data_path = str(os.path.join(base_path, 'process_data.py'))

run_cmd('python ' + process_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --series_number ' + series_number_for_roi_selection
                + ' --run_gui ' + run_gui
                + ' --fano_regularizer ' + fano_regularizer
                + ' --attach_metadata ' + attach_metadata
                + ' --roi_set_name ' + roi_set_name
                + ' --response_set_name_prefix ' + response_set_name
                + ' --flymax_movies_root ' + flymax_movies_root)

# %%
