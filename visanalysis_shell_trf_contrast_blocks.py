"""
visanalysis shell script for running analyze_data_trf_contrast_blocks.py:
TEMPORAL filters under switching contrast (ContrastBlockUniformNoise), one set
per roi per channel, on flies whose rois already exist.

These recordings (Mi1 split, JS140_x_JS252 = gflamp2 + geco; ch1 = calcium,
ch2 = cAMP) were processed in July 2025 with process_data/select_rois, so there
is no processing block here -- only the filter analysis. Each fly has:
    series 1  search stimulus (MovingPatch)
    series 2  contrast blocks, FROZEN noise (start_seed 0, 60 epochs)
    series 3  contrast blocks, random seeds (300 epochs, ~25 min)
The script tells them apart by seed, not by number: series 3 gives the
filters, series 2 the noise ceiling and an out-of-series test of them.

`experiment_file_directories` is a LIST; every block loops over it.

Run top to bottom, or re-run the #%% cells individually.
"""
#%% INITIALIZE ENVIRONMENT

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(
    globals().get('__file__', 'visanalysis_shell_trf_contrast_blocks.py'))))
from visanalysis.util.cellrun import run   # prints child output into the cell


#%% INITIALIZE ARGUMENTS

base_path = 'C:/Users/jcsimon/Documents/GitHub/visanalysis'

experiment_file_directories = [
    # fly_011 excluded (2026-09-23): no final ROI selection, 11/57 ROIs pass, noisy
    # 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS140_x_JS252/fly_011',
    'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS140_x_JS252/fly_012',
    'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS140_x_JS252/fly_013',
]

rig = 'Bruker'
tag = 'final'                # reads fly_final.hdf5 where it exists, else fly.hdf5
                             # (fly_011 never had final rois picked; the null
                             # gate does the selecting there)
show_figs = 'False'
save_figs = 'True'

filter_length = '2.0'        # s. Mi1 calcium filters are over by ~0.5 s
steady_start = ''            # s since switch; '' = filter_length
gain_bins = '0,0.1,0.2,0.3,0.5,0.75,1,1.5,2,3,4,5.1'
null_shifts = '20'
lag_dt = '0.01'              # s. Sub-volume lag grid: samples fall at random
                             # phases of the 30 Hz stimulus, so lags finer than
                             # the 72 ms volume are resolvable ('0' = volume)
sample_window = 'frame'      # each plane integrates one frame period (~5.5 ms)
channel_names = '1:calcium (GECO),2:cAMP (gFlamp2)'   # JS140_x_JS252
sta_bin = '0.02'             # s, jitter-binned switch-triggered averages
channels = '1,2'

# '' writes into each fly's hdf5 (group TRF_CONTRAST_BLOCKS_<tag>) and its
# <tag>_roi_figs folder. A path writes there instead, one subfolder per fly --
# use it when the data folder is read-only to the account running this.
output_root = 'C:/Users/Public/flymax_analysis/contrast_blocks'


#%% ANALYZE: filters under switching contrast

from visanalysis.util.cellrun import run

for experiment_file_directory in experiment_file_directories:
    print('\n=== analyze_data_trf_contrast_blocks: ' + experiment_file_directory
          + ' ===')
    cmd = ['--experiment_file_directory', experiment_file_directory,
           '--rig', rig,
           '--show_figs', show_figs,
           '--save_figs', save_figs,
           '--tag', tag,
           '--filter_length', filter_length,
           '--gain_bins', gain_bins,
           '--null_shifts', null_shifts,
           '--channels', channels,
           '--lag_dt', lag_dt,
           '--sample_window', sample_window,
           '--sta_bin', sta_bin,
           '--channel_names', channel_names]
    if steady_start:
        cmd += ['--steady_start', steady_start]
    if output_root:
        fly = os.path.basename(os.path.normpath(experiment_file_directory))
        geno = os.path.basename(os.path.dirname(
            os.path.normpath(experiment_file_directory)))
        cmd += ['--output_directory', os.path.join(output_root, geno, fly)]
    run(os.path.join(base_path, 'analyze_data_trf_contrast_blocks.py'), *cmd)
