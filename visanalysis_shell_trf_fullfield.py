"""
visanalysis shell script for running:
1. process_data.py (series 1: identify rois from search/flash stim)
2. analyze_data.py (raw)
3. select_rois.py
4. analyze_data.py (final)
5. analyze_data_trf_fullfield.py (final) - TEMPORAL filters from the full-field
   ternary noise series, one per roi, using final rois

Same block structure as visanalysis_shell_strf_octave.py, with one difference:
`experiment_file_directories` is a LIST. Every block loops over it, so several
flies run from one execution -- which is the normal case for this stimulus,
since a single 300 s trial per fly is enough for a temporal filter and the
interesting comparison is across flies.

Run it top to bottom the first time, then re-run individual blocks (#%% cells)
while tuning. The GUI blocks (process_data with run_gui, select_rois) stop for
each fly in turn.

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
#%% INITIALIZE ENVIRONMENT

import os
import sys

# Import, do not define: a helper defined in this cell disappears whenever the
# file is edited mid-session, and every later block then fails with
# "NameError: name 'run' is not defined". Each block imports it too, so blocks
# depend only on the ARGUMENTS cell.
sys.path.insert(0, os.path.dirname(os.path.abspath(
    globals().get('__file__', 'visanalysis_shell_trf_fullfield.py'))))
from visanalysis.util.cellrun import run   # prints child output into the cell


#%% INITIALIZE ARGUMENTS

# all scripts
base_path = 'C:/Users/jcsimon/Documents/GitHub/visanalysis'

# one entry per fly folder containing a fly.hdf5. Add or comment out lines here;
# nothing below needs editing to change which flies run.
experiment_file_directories = [
    'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS257/fly_014',
]

rig = 'Bruker' #string "Bruker" or "AODscope"

roi_set_name = 'roi_set_name' # name of roi group to be analyzed (default 'roi_set_name')
response_set_name = 'mask' # name of response group to be analyzed (default 'mask')

# process_data
series_number_for_roi_selection = '2' #string, first series (search/flash stim) is used to identify rois
run_gui = 'True' #string "True" or "False", default = "False"
fano_regularizer = '300' #string, k in the GUI Fano image var / (mean + k), in image intensity units
attach_metadata = 'True' #string "True" or "False", default = "False"
flymax_movies_root = 'C:/Users/jcsimon/Documents/Stanford/Data/flymax_movies' #root dir to search for noise stimulus .bin/_params.json files
                                                                              #unused by the full-field stimulus, which has no movie file,
                                                                              #but process_data still takes it for other series

# analyze_data
show_figs = 'False' #string "True" or "False", default = "False"
save_figs = 'True' #string "True" or "False", default = "False"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

# select_rois
save_hdf5 = 'True' #string "True" or "False", default = "False"

# analyze_data_trf_fullfield
filter_length = '10' #string, filter length in seconds. Must cover the slowest kernel you
                     #expect -- cAMP indicators still carry power past 3 s. Too short shows
                     #up as filters that have not returned to zero, which the script counts
block_s = '60' #string, cross-validation block length in seconds. Folds are contiguous blocks
               #rather than epochs, because the protocol ships ONE 300 s trial and a fold per
               #epoch would leave nothing to hold out. Shorter gives more folds and less
               #training data each; 60 s is 5 folds per trial
ridge = '' #string, fix the ridge instead of cross-validating it. Leave '' -- the right value
           #depends on the hold time and on the indicator, so a number that suited one
           #recording will not suit the next. Set it only to reproduce an old analysis
band = '' #string, fit against ONE band of a multi-band stimulus instead of the displayed
          #composite. Leave '' for the normal single-band recording. Diagnostic only: the
          #composite is what the cell actually saw
bleach = 'auto' #string "auto", "on" or "off". 'auto' corrects when the epoch is long enough
                #to separate bleaching from signal and skips it otherwise
null_shifts = '20' #string, circular-shift null draws; '0' disables. Costs (n+1)x the fit.
                   #This is what decides which rois are reported as driven at all.
                   #Below 20, FDR cannot pass a single responsive roi at q=0.05
                   #(min_detectable_k ~= 1/(0.05*n_shifts)), so both gates fall
                   #back to a percentile rule estimated from few draws.


#%% PROCESS_DATA
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run


process_data_path = str(os.path.join(base_path, 'process_data.py'))

for experiment_file_directory in experiment_file_directories:
    print('\n=== process_data: ' + experiment_file_directory + ' ===')
    run(process_data_path,
        '--experiment_file_directory', experiment_file_directory,
        '--rig', rig,
        '--series_number', series_number_for_roi_selection,
        '--run_gui', run_gui,
        '--fano_regularizer', fano_regularizer,
        '--attach_metadata', attach_metadata,
        '--roi_set_name', roi_set_name,
        '--response_set_name_prefix', response_set_name,
        '--flymax_movies_root', flymax_movies_root)


#%% ANALYZE_DATA RAW
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run


tag = 'raw' #string "raw" or "final"

analyze_data_path = str(os.path.join(base_path, 'analyze_data.py'))

for experiment_file_directory in experiment_file_directories:
    print('\n=== analyze_data raw: ' + experiment_file_directory + ' ===')
    run(analyze_data_path,
        '--experiment_file_directory', experiment_file_directory,
        '--rig', rig,
        '--show_figs', show_figs,
        '--save_figs', save_figs,
        '--tag', tag,
        '--dff', dff)


#%% SELECT_FINAL_ROIS
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run
input_tag = ''
output_tag = 'final'
select_rois_path = str(os.path.join(base_path, 'select_rois.py'))

for experiment_file_directory in experiment_file_directories:
    print('\n=== select_rois: ' + experiment_file_directory + ' ===')
    run(select_rois_path,
        '--experiment_file_directory', experiment_file_directory,
        '--rig', rig,
        '--save', save_hdf5,
        '--input_tag', input_tag,
        '--output_tag', output_tag)


#%% ANALYZE_DATA FINAL
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run


tag = 'final' #string "raw" or "final"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

analyze_data_path = str(os.path.join(base_path, 'analyze_data.py'))

for experiment_file_directory in experiment_file_directories:
    print('\n=== analyze_data final: ' + experiment_file_directory + ' ===')
    run(analyze_data_path,
        '--experiment_file_directory', experiment_file_directory,
        '--rig', rig,
        '--show_figs', show_figs,
        '--save_figs', save_figs,
        '--tag', tag,
        '--dff', dff)


#%% ANALYZE_DATA_TRF_FULLFIELD FINAL
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run
# uses the final rois selected above; auto-detects the full-field ternary series
# within each fly by protocol_ID, so no series number is hard-coded here.
#
# As in the octave path there is no movie to locate: the stimulus is regenerated
# from the epoch metadata, so this block needs nothing on disk beyond the hdf5.
#
# There is deliberately no --dff here. The filter is fitted against dF/F taken
# from the per-roi MEDIAN, not from the pre period: two seconds of baseline in
# front of a 300 s epoch is a noisy divisor, and the median over the epoch
# measures the same quantity ~150x better.
#
# There is also no --deconvolve. The octave path needs one because an STA
# returns the filter convolved with the stimulus autocorrelation; this is a
# regression, so that is divided out as part of the fit and the peak lag is
# already the peak lag.

tag = 'final' #string "raw" or "final"

analyze_data_trf_path = str(os.path.join(base_path, 'analyze_data_trf_fullfield.py'))

for experiment_file_directory in experiment_file_directories:
    print('\n=== analyze_data_trf_fullfield: ' + experiment_file_directory + ' ===')
    args = ['--experiment_file_directory', experiment_file_directory,
            '--rig', rig,
            '--show_figs', show_figs,
            '--save_figs', save_figs,
            '--tag', tag,
            '--filter_length', filter_length,
            '--block_s', block_s,
            '--bleach', bleach,
            '--null_shifts', null_shifts]
    if ridge != '':
        args += ['--ridge', ridge]
    if band != '':
        args += ['--band', band]

    run(analyze_data_trf_path, *args)
# %%
