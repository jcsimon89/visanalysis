"""
visanalysis shell script for running analyze_data_trf_contrast_switch.py:
contrast adaptation from ContrastSwitchTernaryNoise, on flies whose ROIs
already exist (process_data / select_rois). The contrast-switch series is found
by protocol_ID, not by number.

Tm9 (JS256 x JS261 control, JS256 x JS257 PDE): ch1 = calcium (GECO),
ch2 = cAMP. The gate runs on the first channel listed, and the ROIs it keeps
are used for every channel.

`experiment_file_directories` is a LIST; every block loops over it.

Run top to bottom, or re-run the #%% cells individually.
"""
#%% INITIALIZE ENVIRONMENT

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(
    globals().get('__file__', 'visanalysis_shell_trf_contrast_switch.py'))))
from visanalysis.util.cellrun import run   # prints child output into the cell


#%% INITIALIZE ARGUMENTS

base_path = 'C:/Users/jcsimon/Documents/GitHub/visanalysis'

experiment_file_directories = [
    'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS261/fly_011',
    # recorded as JS256_x_JS261 but is JS256_x_JS257 (PDE); folder not yet renamed
    'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS261/fly_012',
]

rig = 'Bruker'
tag = 'raw'                  # reads fly_<tag>.hdf5 where it exists, else fly.hdf5
show_figs = 'False'
save_figs = 'True'

channels = '1:calcium,2:cAMP'          # first = primary: the gate runs on it
filter_length = '1:3.0,2:5.0'          # s
lag_dt = '1:0.02,2:0.1'                # s; finer than the volume, from the plane timing
steady_start = ''                      # s after a switch; '' = the filter length
first_dwell = '2'
gate = 'high'                          # or a path to a .npy of ROI indices
gate_q = '0.05'
null_rolls = '14'
n_boot = '1000'
crossfit = 'True'
export_models = 'False'   # True: write each channel's steady-state LN pair to the model library
model_library = ''        # '' = $VISANALYSIS_MODEL_LIBRARY or ~/.visanalysis/config.yaml

# '' writes into <fly folder>/<tag>_contrast_switch/. A path writes there
# instead, one subfolder per fly -- use it when the data folder is read-only to
# the account running this.
output_root = 'C:/Users/Public/flymax_analysis/contrast_switch'


#%% ANALYZE

for experiment_file_directory in experiment_file_directories:
    print('\n=== analyze_data_trf_contrast_switch: ' + experiment_file_directory + ' ===')
    cmd = ['--experiment_file_directory', experiment_file_directory,
           '--rig', rig,
           '--show_figs', show_figs,
           '--save_figs', save_figs,
           '--tag', tag,
           '--channels', channels,
           '--filter_length', filter_length,
           '--lag_dt', lag_dt,
           '--first_dwell', first_dwell,
           '--gate', gate,
           '--gate_q', gate_q,
           '--null_rolls', null_rolls,
           '--n_boot', n_boot,
           '--crossfit', crossfit,
           '--export_models', export_models]
    if model_library:
        cmd += ['--model_library', model_library]
    if steady_start:
        cmd += ['--steady_start', steady_start]
    if output_root:
        fly = os.path.basename(os.path.normpath(experiment_file_directory))
        geno = os.path.basename(os.path.dirname(os.path.normpath(experiment_file_directory)))
        cmd += ['--output_directory', os.path.join(output_root, geno, fly)]
    run(os.path.join(base_path, 'analyze_data_trf_contrast_switch.py'), *cmd)
