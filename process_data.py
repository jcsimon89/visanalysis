"""
process data:
1. attach metadata
2. draw roi masks with reduced gui
3. extract roi responses and save to hdf5

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""

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
from visanalysis.util import plot_tools
from visanalysis.util import h5io
from visanalysis.util.noise_stim import (
    is_flymax_movie, get_epoch_bin_filename, find_bin_files, load_stim_params_json,
    resolve_attrs_to_write, unique_values,
)

# call structure: python process_data.py --experiment_file_directory "path" --rig "rigID" --series_number "series_number" --run_gui "True/False" --attach_metadata "True/False"

# experiment_file_directory: string that contains the path to the folder that has the .hdf5 file
# rig: rigID (Bruker or AODscope)

#NOTE: all file names and paths are relative to the experiment_file_directory which is the main fly folder ie /fly_001


if __name__ == '__main__':
    ## parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment_file_directory", nargs="?", help="Folder pointing to hdf5")
    parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
    parser.add_argument("--series_number", nargs="?", help="number as string for series to use for roi selection")
    parser.add_argument("--run_gui", nargs="?", help="True/False")
    parser.add_argument("--attach_metadata", nargs="?", help="True/False")
    parser.add_argument("--roi_set_name", nargs="?", default='roi_set_name', help="name of roi set for analysis")
    parser.add_argument("--response_set_name_prefix", nargs="?", default='mask', help="name of response set for analysis")
    parser.add_argument("--flymax_movies_root", nargs="?", default=None,
                        help="Root directory to search for noise stimulus .bin/_params.json files "
                             "(only required if the fly has a noizone series)")
    args = parser.parse_args()

    roi_set_name = args.roi_set_name
    response_set_name_prefix = args.response_set_name_prefix
    experiment_file_directory = args.experiment_file_directory
    rig = args.rig
    series_number = args.series_number
    flymax_movies_root = args.flymax_movies_root
    
    if args.run_gui == 'True':
        run_gui = True
    elif args.run_gui == 'False':
        run_gui = False
    else:
        print('not able to interperet run_gui flag, must be "True" or "False"')
    
    if args.attach_metadata == 'True':
        attach_metadata = True
    elif args.attach_metadata == 'False':
        attach_metadata = False
    else:
        print('not able to interperet attach_metadata flag, must be "True" or "False"')


    # hardcoded names
    experiment_file_name = 'fly.hdf5' #name of hdf5 file
    json_file_name = 'fly.json' #name of json file

    # extract info from fly.json 

    with open(pathlib.Path(experiment_file_directory, json_file_name), 'r') as file:
        fly_json = json.load(file)

    struct_channel = [fly_json['structural_channel']] # can only be one structural channel
    struct_channel_num = [struct_channel[0].split('_')[-1]] #datatype=list of strings
    print('struct_channel = ' + str(struct_channel))
    print('struct_channel_num = ' + str(struct_channel_num))

    func_channels = fly_json['functional_channel'].replace("[","").replace("]","").replace("'","").split(",") # can be many functional channels - weird format in snake_brainsss json imports as one big string, cleaning up and converting to list of strings #datatype=list of strings
    func_channels_num = [func_channels[i].split('_')[-1] for i in range(len(func_channels))] #datatype = list of strings, but individual channel numbers will be converted to int before using methods
    print('func_channels = ' + str(func_channels))
    print('func_channel_num = ' + str(func_channels_num))

    #TODO: other data needed from fly_json?
    DO_BG_CH1 = fly_json.get('background_subtraction_ch1', False)
    DO_BG_CH2 = fly_json.get('background_subtraction_ch2', False)

    # derive image_file_path
    # assumes folder structure from snake_eyesss
    _struct_ch = struct_channel_num[0]
    _do_bg_struct = (DO_BG_CH1 if _struct_ch == '1' else DO_BG_CH2 if _struct_ch == '2' else False)
    if _do_bg_struct:
        image_file_name = 'channel_' + _struct_ch + '_moco_bg_func.nii'
    else:
        image_file_name = 'channel_' + _struct_ch + '_moco_func.nii'

    image_relative_directory = 'func' + str(int(series_number)-1) + '/moco' #folder where .nii is, assumes func_ folder counting starts from 0 which series counter starts from 1
    image_file_directory = os.path.join(experiment_file_directory, image_relative_directory)
    image_file_path = os.path.join(image_file_directory, image_file_name)

    print('Image file location for roi selection: ' + str(image_file_path))

    # file_path (complete path to .hdf5 file)
    experiment_file_path = os.path.join(experiment_file_directory, experiment_file_name)

    #TODO extract subject number from hdf5 or folder name? not sure if needed

    print('experiment_file_path: ' + repr(os.path.join(experiment_file_directory, experiment_file_name)))

    ## determine correct visanalysis plugin

    if rig == 'Bruker':
        from visanalysis.plugin import bruker
        plug = bruker.BrukerPlugin()
        print('****Bruker plugin****')
    elif rig == 'AODscope':
        from visanalysis.plugin import aodscope
        plug = aodscope.AodScopePlugin()
        print('****AODscope plugin****')
    else:
        plug = base_plugin.BasePlugin()
        print('****Unrecognized plugin name****')

    # attach metadata to hdf5
    print('attach_metadata: ' + str(attach_metadata))
    if attach_metadata:
        plug.attachData(experiment_file_name, experiment_file_directory)

        print('Attached metadata to {}'.format(experiment_file_name))

    # ── resolve + cache movie-stimulus parameters (any flymax movie-playback ──
    # series -- noizone, the Edge/search stimulus, or anything else built the
    # same way) ─────────────────────────────────────────────────────────────
    # For each such series, look up every epoch's bin filename, resolve all of
    # them under flymax_movies_root in one pass, load each unique bin's
    # _params.json sidecar once, and write bin_path + all json params onto
    # each epoch as hdf5 attributes -- so downstream analysis reads them
    # straight from the hdf5 instead of touching the filesystem. Also writes,
    # at the series (run) level, the set of unique values seen for each field
    # across that series' epochs -- mirroring what native stimpack itself
    # writes for a live-rendered protocol (e.g. run_parameters['intensity']
    # listing every value used across the series). Before writing anything,
    # every value is checked against what's already there (native stimpack
    # fields, or a previous run of this script) -- same value: skip; missing:
    # write; different value: hard error, since that means a same-named field
    # means something different natively and would otherwise be silently
    # clobbered (see resolve_attrs_to_write).
    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path)))
    print('series_num = ' + str(series_num))

    movie_epochs_by_series = {}  # series_int -> list of (epoch_ind, bin_filename)
    run_params_by_series = {}    # series_int -> run_parameters (for the collision check)
    needed_filenames = set()

    for current_series in series_num:
        current_series_int = int(current_series)
        plug.updateImagingDataObject(experiment_file_directory, experiment_file_name, current_series_int)
        movie_ID = plug.ImagingDataObject
        run_params = movie_ID.getRunParameters()
        if not is_flymax_movie(run_params):
            continue
        run_params_by_series[current_series_int] = run_params
        epoch_list = []
        for epoch_ind, epoch_params in enumerate(movie_ID.getEpochParameters()):
            bin_filename = get_epoch_bin_filename(epoch_params)
            epoch_list.append((epoch_ind, bin_filename))
            needed_filenames.add(bin_filename)
        movie_epochs_by_series[current_series_int] = epoch_list

    if movie_epochs_by_series:
        if not flymax_movies_root:
            raise ValueError(
                'Found flymax movie-playback series but --flymax_movies_root was not '
                'provided. Pass the root directory containing the movie .bin/_params.json files.'
            )
        print('Resolving {} unique movie bin file(s) under {}...'.format(
            len(needed_filenames), flymax_movies_root))
        bin_paths = find_bin_files(flymax_movies_root, needed_filenames)

        params_cache = {}  # bin_path -> params dict (loaded once per unique file)
        for bin_path in set(bin_paths.values()):
            params_json = load_stim_params_json(bin_path)
            if params_json is None:
                raise FileNotFoundError('No _params.json sidecar found next to {}'.format(bin_path))
            params_cache[bin_path] = params_json

        def _cache_movie_attrs(target_file_name, target_file_path):
            """
            Resolve + write the movie attrs against target_file_path specifically
            (its own current attrs, via its own ImagingDataObject) rather than
            reusing the raw file's collision check -- so this can be called again
            for fly_final.hdf5 (once select_rois.py has created it) to keep it in
            sync without ever touching its ROI/STRF datasets.
            """
            n_epochs = 0
            for current_series_int, epoch_list in movie_epochs_by_series.items():
                resolved_epoch_attrs = []  # each epoch's full resolved attrs, for the run-level uniques below

                plug.updateImagingDataObject(experiment_file_directory, target_file_name, current_series_int)
                target_epoch_params_list = plug.ImagingDataObject.getEpochParameters()
                target_run_params = plug.ImagingDataObject.getRunParameters()

                for epoch_ind, bin_filename in epoch_list:
                    bin_path = bin_paths[bin_filename]
                    new_attrs = {'bin_path': bin_path, **params_cache[bin_path]}
                    to_write = resolve_attrs_to_write(
                        target_epoch_params_list[epoch_ind], new_attrs,
                        '{} series {} epoch {}'.format(target_file_name, current_series_int, epoch_ind),
                    )
                    if to_write:
                        h5io.updateEpochAttributes(target_file_path, current_series_int, epoch_ind, to_write)
                    resolved_epoch_attrs.append(new_attrs)
                    n_epochs += 1

                # run-level: unique values per field across all epochs in this series
                all_keys = set()
                for attrs in resolved_epoch_attrs:
                    all_keys.update(attrs.keys())
                run_level_attrs = {}
                for key in all_keys:
                    values = unique_values([a[key] for a in resolved_epoch_attrs if key in a])
                    run_level_attrs[key] = values[0] if len(values) == 1 else np.array(values)

                to_write_run = resolve_attrs_to_write(
                    target_run_params, run_level_attrs,
                    '{} series {} (run-level)'.format(target_file_name, current_series_int),
                )
                if to_write_run:
                    h5io.updateSeriesAttributes(target_file_path, current_series_int, to_write_run)
            return n_epochs

        n_epochs_total = _cache_movie_attrs(experiment_file_name, experiment_file_path)
        print('Cached movie parameters for {} epoch(s) across {} movie-playback series.'.format(
            n_epochs_total, len(movie_epochs_by_series)))

        # Mirror the same attrs into fly_final.hdf5 if select_rois.py has already
        # created it (a prior pipeline run) -- keeps it in sync with fly.hdf5
        # without ever re-copying/touching its ROI or STRF datasets. See
        # select_rois.py: once fly_final.hdf5 exists, it only ever patches ROI
        # datasets in place and never re-syncs attrs from the raw file itself.
        final_file_name = 'fly_final.hdf5'
        final_file_path = os.path.join(experiment_file_directory, final_file_name)
        if os.path.exists(final_file_path):
            print('{} already exists -- mirroring cached movie parameters into it too...'.format(final_file_name))
            _cache_movie_attrs(final_file_name, final_file_path)

    ##draw roi masks using reduced GUI
    print('run_gui: ' + str(run_gui))
    if run_gui:

        gui_path = str(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gui", "DataGUI_prog.py"))

        os.system('python ' + gui_path
                + ' --experiment_file_directory ' + experiment_file_directory
                + ' --experiment_file_name ' + experiment_file_name
                + ' --experiment_file_path ' + experiment_file_path
                + ' --rig ' + rig
                + ' --series_number ' + series_number
                + ' --image_file_path ' + image_file_path)


    ID = imaging_data.ImagingDataObject(experiment_file_path,
                                        series_number,
                                        quiet=False)

    # Retrieve saved mask region responses from data file

    roi_mask_bool = ID.getRoiMasks(roi_set_name)['roi_mask'] #shape: roi_index, x, y ,(z)
    roi_image = ID.getRoiMasks(roi_set_name)['roi_image'] #shape: x,y,(z)

    print('dimensions of roi_image: ' + str(np.shape(roi_image)))
    print('dimensions of roi_mask_bool: ' + str(np.shape(roi_mask_bool)))

    # figure out if data is volume or slice
    num_spacial_dim = len(roi_image.shape)

    # add func_spatial_dim and get fly metadata for later
    fly_metadata = ID.getSubjectMetadata()
    print('fly_metadata: ' + repr(fly_metadata))
   
    # convert roi mask to unique numbers for each separate roi (background = 0, first roi =1, second roi =2 etc.)
    # this is how masks need to be formatted for plug.saveRegionResponsesFromMask()
    roi_mask = np.zeros(roi_mask_bool.shape[1:])
    print('dimensions of func data: ' + str(np.shape(roi_mask))) # shape:x,y,(z)

    if num_spacial_dim == 3: #data is a volume
        for roi_ind in range(roi_mask_bool.shape[0]):
            for i in range(roi_mask_bool.shape[1]):
                for j in range(roi_mask_bool.shape[2]):
                    for k in range(roi_mask_bool.shape[3]):
                        if roi_mask_bool[roi_ind,i,j,k] == True:
                            roi_mask[i,j,k]=roi_ind+1 #since roi_ind starts at 0 and we want to label the first roi with value=1

    elif num_spacial_dim == 2: #data is a slice
        for roi_ind in range(roi_mask_bool.shape[0]):
            for i in range(roi_mask_bool.shape[1]):
                for j in range(roi_mask_bool.shape[2]):
                    if roi_mask_bool[roi_ind,i,j] == True:
                            roi_mask[i,j]=roi_ind+1 #since roi_ind starts at 0 and we want to label the first roi with value=1
    else:
        print('data does not appear to be a volume or a slice, num_spatial_dim = ' + str(num_spacial_dim))

    n_roi = len(np.unique(roi_mask))-1 # numpy integer
    print('number of rois in roi_mask: ' + str(n_roi))


    ## start response extraction
    # (series_num already computed above, before ROI selection)

    for series_ind, current_series in enumerate(series_num): #loop through all series
        current_series = int(current_series) # methods expect series number to be datatype int
        sn = 'sn' + str(current_series)

        #update imaging object with current series number

        plug.updateImagingDataObject(experiment_file_directory, #NOTE: does this set current_series in imaging data object? I think so (instantiates self.ImagingObject with current_series)
                                    experiment_file_name,
                                    current_series)
        
        ID = plug.ImagingDataObject
    
        for current_channel in func_channels_num: #loop through channels
            current_channel = int(current_channel) # methods expect channel number to be datatype int
            ch = 'ch' + str(current_channel)
            response_set_name = response_set_name_prefix + '_' + ch

            #derive image file name — moco aligns all series to func0's template directly
            _do_bg = (DO_BG_CH1 if current_channel == 1 else DO_BG_CH2 if current_channel == 2 else False)
            if _do_bg:
                image_file_name = 'channel_' + str(current_channel) + '_moco_bg_func.nii'
            else:
                image_file_name = 'channel_' + str(current_channel) + '_moco_func.nii'


            image_relative_directory = 'func' + str(int(current_series)-1) + '/moco' #folder where .nii is, assumes func_ folder counting starts from 0 which series counter starts from 1
            image_file_directory = os.path.join(experiment_file_directory, image_relative_directory)
            image_file_path = os.path.join(image_file_directory, image_file_name)
            
            #associate image data

            plug.updateImageSeries(data_directory=image_file_directory, #NOTE: doesnt do anything with current_series!
                                    image_file_name=image_file_name,
                                    series_number=current_series,
                                    channel=current_channel)

            # save mask and responses to hdf5 (roi_prefix = aligned)

            plug.saveRegionResponsesFromMask(file_path=experiment_file_path,
                                                series_number=current_series,
                                                response_set_name=response_set_name,
                                                mask=roi_mask,
                                                include_zero=False)
            
            