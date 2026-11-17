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
experiment_file_head = 'C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/'
experiment_file_tail = ['JS140_x_JS261/fly_023']


# EVERY FLY THE STRF CELL RUNS. One entry is the single-fly behaviour; add paths
# to batch. The stimulus is regenerated from each fly's own epoch metadata, so
# nothing here has to match between them, and one fly failing does not lose the
# rest -- the exit code is collected and a summary names the failures.
#
# A fly whose saved results already match the settings below is REUSED: the
# strfs are read back and only the figures are redrawn, a few minutes rather
# than a recompute. Only filter_length, dff, deconvolve and lag_dt affect that
# decision; every display option is free to change. --force_reanalysis overrules.
# head + each tail. A list cannot be concatenated to a string, which is what
# `experiment_file_head + experiment_file_tail` was trying to do -- the join has
# to happen once per entry.
fly_dirs = [experiment_file_head + t for t in experiment_file_tail]

# The other cells still take ONE fly, so they default to the first of the list.
# Point it somewhere else if you want to run those against a different fly.
experiment_file_directory = fly_dirs[0]
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

# analyze_data
show_figs = 'False' #string "True" or "False", default = "False"
save_figs = 'True' #string "True" or "False", default = "False"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

# select_final_rois
save_hdf5 = 'True' #string "True" or "False", default = "False"

# analyze_data_strf_one_d
filter_length = '10' #string. ONE number for every channel, or PER CHANNEL as
                     #'ch1:3,ch2:25' -- a slow indicator needs a long filter while a fast
                     #one only pays for it, in compute and in a lag axis that is mostly
                     #empty. Must cover the slowest kernel you
                    #expect, on top of the stimulus's own hold: the 0.45 s default single hold
                    #has decayed to under 1% of its autocorrelation by about 2.1 s, so 10 s is
                    #ample headroom for a slow indicator. Lower it for a fast one -- the lag
                    #axis is sampled at the raw display rate, so length costs compute directly
monitor_hz = '' #string, display rate override in Hz, leave '' to use the ImagingDataObject default (120 Hz)
lam = '' #string, deconvolution ridge. Leave '' to cross-validate it across trial halves,
         #which is the intended path -- a fixed value that suits calcium will not suit cAMP
lam_scale = '1.0' #string, multiplies the cross-validated lambda. >1 is smoother, <1 sharper
                  #and noisier. Use this to tune rather than pinning lam outright
radial_symmetry = 'True' #string "True" or "False". Recovers a radial profile h(r)
                         #from the registered means by forward-modelling the strip
                         #projections and solving ONE h against every family at once,
                         #then renders it as a 2D map. Needs two families about
                         #different poles, which a swept pole gives you.
                         #
                         #Read the residual it prints before believing the 2D picture:
                         #with one family the assumption cannot even be tested, and
                         #with two the residual is what falsifies it. Above ~0.15 the
                         #field is not radially symmetric and the render is
                         #illustrating the assumption rather than the measurement.
mean_span_deg = '' #string, half-width in degrees of the registered-mean
                   #figures. '' uses 75. The registered grid runs as far as
                   #the widest roi reached, which on a wedges family is most
                   #of a hemisphere; clamping is what makes the mean filter
                   #readable rather than a few pixels wide
n_display_lags = '' #string, how many lags the 2D lag grid shows, spread over the
                    #filter support. '' uses 12. One frame at the peak lag is not
                    #enough when center and surround peak at different times
figure_rois = '' #string 'kept' (default), 'none', 'all', or a number N for the N strongest
                 #kept rois. '' uses 'kept'. 'all' on a 289-roi fly with two
                 #geometries and two channels is 2312 files, most of them from a
                 #channel where nothing passed
lag_dt = '' #string, lag axis step in SECONDS. '' uses 1/60 s, half the display rate.
            #One lag per display frame is ~14 bins per response sample at this rig's
            #imaging rate, which costs compute without resolving anything

deconvolve = 'False' #string "True" or "False". OFF by default, and that is the right default:
                     #the hold-smear leaves the POSITIONAL profile exact for a separable filter,
                     #so the bar/wedge centroid and width need no correction. It biases TIME
                     #only -- peak lag late by roughly hold_tau -- and that bias is nearly
                     #constant across rois, so latency DIFFERENCES between them are already
                     #unbiased. Set "True" only if you need absolute latency, then read the
                     #printed decision: the reliability guard refuses it where it would cost
                     #more than it buys


# #%% PROCESS_DATA
# import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
# from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell


# process_data_path = str(os.path.join(base_path, 'process_data.py'))

# run_cmd('python ' + process_data_path
#                 + ' --experiment_file_directory ' + experiment_file_directory
#                 + ' --rig ' + rig
#                 + ' --series_number ' + series_number_for_roi_selection
#                 + ' --run_gui ' + run_gui
#                 + ' --fano_regularizer ' + fano_regularizer
#                 + ' --attach_metadata ' + attach_metadata
#                 + ' --roi_set_name ' + roi_set_name
#                 + ' --response_set_name_prefix ' + response_set_name
#                 + ' --flymax_movies_root ' + flymax_movies_root)


#%% ANALYZE_DATA RAW
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell


tag = 'raw' #string "raw" or "final"

analyze_data_path = str(os.path.join(base_path, 'analyze_data.py'))

run_cmd('python ' + analyze_data_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff)


# #%% SELECT_FINAL_ROIS
# import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
# from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell
# input_tag = ''
# output_tag = 'final'
# select_rois_path = str(os.path.join(base_path, 'select_rois.py'))

# run_cmd('python ' + select_rois_path
#                 + ' --experiment_file_directory ' + experiment_file_directory
#                 + ' --rig ' + rig
#                 + ' --save ' + save_hdf5
#                 + ' --input_tag ' + input_tag
#                 + ' --output_tag ' + output_tag)


# #%% ANALYZE_DATA FINAL
# import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
# from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell


# tag = 'final' #string "raw" or "final"
# dff = 'pre' #string "pre", "mean", or "none", default = "pre"

# analyze_data_path = str(os.path.join(base_path, 'analyze_data.py'))

# run_cmd('python ' + analyze_data_path
#                 + ' --experiment_file_directory ' + experiment_file_directory
#                 + ' --rig ' + rig
#                 + ' --show_figs ' + show_figs
#                 + ' --save_figs ' + save_figs
#                 + ' --tag ' + tag
#                 + ' --dff ' + dff)


#%% ANALYZE_DATA_STRF_ONE_D FINAL
import sys; sys.path.insert(0, base_path) if base_path not in sys.path else None
from visanalysis.util.cellrun import run_cmd   # prints output into THIS cell
# uses the final rois selected above; auto-detects the 1D ternary series within
# the fly by protocol_ID, so the series number is not hard-coded here.
#
# Like the octave path there is no movie to locate: the stimulus is regenerated
# from the epoch metadata, so this block needs nothing on disk beyond the hdf5.

tag = 'raw' #string "raw" or "final"
dff = 'pre' #string "pre", "mean", or "none", default = "pre"

majority_sign_only = 'True' #string "True" or "False". Adds a mean over only the
                             #rois carrying the channel's majority polarity. Fair when
                             #the sign is known from the cell type (Mi1 calcium is ON)
                             #rather than read off these data. NOT sign-alignment --
                             #dissenters are dropped, never flipped.
edge_margins = '' #string, e.g. '0,10,20'. One extra registered mean per margin,
                  #from the rois at least that far inside the OBSERVED extent of
                  #the population's rf centers. '' is off.
                  #
                  #Read them as a series rather than picking one. If the kinetics
                  #do not move between margins, edge rois are not affecting them
                  #and the full population can be kept on evidence. A fixed margin
                  #does NOT mean the same thing on every fly: it bites according
                  #to how compact that fly's population is, and 20 deg left 32% of
                  #fly_016 but only 11% of fly_010.

analyze_data_strf_path = str(os.path.join(base_path, 'analyze_data_strf_one_d.py'))

failed = []
for _i, _fly in enumerate(fly_dirs):
    print('')
    print('=' * 78)
    print('[{}/{}] {}'.format(_i + 1, len(fly_dirs), _fly))
    print('=' * 78)

    command = ('python ' + analyze_data_strf_path
                + ' --experiment_file_directory ' + _fly
                + ' --rig ' + rig
                + ' --show_figs ' + show_figs
                + ' --save_figs ' + save_figs
                + ' --tag ' + tag
                + ' --dff ' + dff
                + ' --filter_length "' + filter_length + '"'
                + ' --lam_scale ' + lam_scale)
    if monitor_hz != '':
        command += ' --monitor_hz ' + monitor_hz
    if lam != '':
        command += ' --lam ' + lam
    if radial_symmetry == 'True':
        command += ' --radial_symmetry'
    if mean_span_deg != '':
        command += ' --mean_span_deg ' + mean_span_deg
    if edge_margins != '':
        command += ' --edge_margins "' + edge_margins + '"'
    if majority_sign_only == 'True':
        command += ' --majority_sign_only'
    if n_display_lags != '':
        command += ' --n_display_lags ' + n_display_lags
    if figure_rois != '':
        command += ' --figure_rois ' + figure_rois
    if lag_dt != '':
        command += ' --lag_dt ' + lag_dt
    if deconvolve == 'True':
        command += ' --deconvolve'

    # One bad fly must not lose the rest of the run: a missing hdf5 or a series
    # that is not 1D ternary is reported and stepped over, not raised.
    _code = run_cmd(command)
    if _code:
        failed.append((_fly, _code))

print('')
print('=' * 78)
print('{} of {} flies completed'.format(len(fly_dirs) - len(failed), len(fly_dirs)))
for _fly, _code in failed:
    print('  FAILED (exit {}): {}'.format(_code, _fly))
# %%
