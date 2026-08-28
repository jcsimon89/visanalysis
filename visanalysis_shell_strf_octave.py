"""
visanalysis shell script for running:
1. process_data.py (series 1: identify rois from search/flash stim)
2. analyze_data.py (raw)
3. select_final_rois.py
4. analyze_data.py (final)
5. analyze_data_strf_octave.py (final) - STRF from the multi-octave ternary
   noise series, in the cell basis, using final rois

Same block structure as visanalysis_shell_strf_noizone.py: run it top to bottom
the first time, then re-run individual blocks (#%% cells) while tuning.

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
#%% INITIALIZE ENVIRONMENT

import os

#%% INITIALIZE ARGUMENTS

# all scripts
base_path = 'C:/Users/jcsimon/Documents/GitHub/visanalysis'
experiment_file_directory = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS140_x_JS261/fly_012' #string to folder containing fly.hdf5 file
rig = 'Bruker' #string "Bruker" or "AODscope"

roi_set_name = 'roi_set_name' # name of roi group to be analyzed (default 'roi_set_name')
response_set_name = 'mask' # name of response group to be analyzed (default 'mask')

# process_data
series_number_for_roi_selection = '1' #string, first series (search/flash stim) is used to identify rois
run_gui = 'True' #string "True" or "False", default = "False"
attach_metadata = 'True' #string "True" or "False", default = "False"
flymax_movies_root = 'C:/Users/jcsimon/Documents/Stanford/Data/flymax_movies' #root dir to search for noise stimulus .bin/_params.json files
                                                                              #unused by the octave stimulus, which has no movie file,
                                                                              #but process_data still takes it for other series

# analyze_data
show_figs = 'False' #string "True" or "False", default = "False"
save_figs = 'True' #string "True" or "False", default = "False"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

# select_final_rois
save_hdf5 = 'True' #string "True" or "False", default = "False"

# analyze_data_strf_octave
filter_length = '10' #string, STRF filter length in seconds. Must cover the slowest kernel
                    #you expect; the 1.00 s coarse hold puts real power out past 3 s
monitor_hz = '' #string, display rate override in Hz, leave '' to use the ImagingDataObject default (120 Hz)
lam = '' #string, deconvolution ridge. Leave '' to cross-validate it across trial halves,
         #which is the intended path -- a fixed value that suits calcium will not suit cAMP
lam_scale = '1.0' #string, multiplies the cross-validated lambda. >1 is smoother, <1 sharper
                  #and noisier. Use this to tune rather than pinning lam outright
deconvolve = 'False' #string "True" or "False". OFF by default, and that is the right default:
                     #the hold-smear leaves the SPATIAL profile exact (measured r = 1.000000),
                     #so centroids, RF size and shape need no correction. It biases TIME only --
                     #peak lag late by 0.08-0.27 s, kernel broadened up to 2.2x -- and that bias
                     #is nearly constant within an octave, so latency DIFFERENCES between rois
                     #are already unbiased. Set "True" only if you need absolute latency, then
                     #read the per-octave decision it prints: the guard refuses it where it
                     #would cost more reliability than it buys


#%% PROCESS_DATA


process_data_path = str(os.path.join(base_path, 'process_data.py'))

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

analyze_data_path = str(os.path.join(base_path, 'analyze_data.py'))

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
select_rois_path = str(os.path.join(base_path, 'select_rois.py'))

os.system('python ' + select_rois_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --save ' + save_hdf5
                + ' --input_tag ' + input_tag
                + ' --output_tag ' + output_tag)


#%% ANALYZE_DATA FINAL


tag = 'final' #string "raw" or "final"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

analyze_data_path = str(os.path.join(base_path, 'analyze_data.py'))

os.system('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff)


#%% ANALYZE_DATA_STRF_OCTAVE FINAL
# uses the final rois selected above; auto-detects the octave ternary series within
# the fly by protocol_ID, so the series number is not hard-coded here.
#
# Unlike the noizone path there is no movie to locate: the stimulus is regenerated
# from the epoch metadata, so this block needs nothing on disk beyond the hdf5.

tag = 'final' #string "raw" or "final"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

analyze_data_strf_path = str(os.path.join(base_path, 'analyze_data_strf_octave.py'))

command = ('python ' + analyze_data_strf_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff
                + ' --filter_length ' + filter_length
                + ' --lam_scale ' + lam_scale)
if monitor_hz != '':
    command += ' --monitor_hz ' + monitor_hz
if lam != '':
    command += ' --lam ' + lam
if deconvolve == 'True':
    command += ' --deconvolve'

os.system(command)
# %%
