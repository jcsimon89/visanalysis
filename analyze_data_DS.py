"""
analyze data:
1. extract metadata from all series
2. extract roi responses from all series, all channels

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
import h5py
import scipy.stats

plt.ioff() #don't show plots unless called for with plt.show() 

# call structure: python analyze_data.py --experiment_file_directory "path" --rig "rigID" --show_figs "True/False" --save_figs "True/False"

# experiment_file_directory: string that contains the path to the folder that has the .hdf5 file

if __name__ == '__main__':
    ## parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment_file_directory", nargs="?", help="Folder pointing to hdf5")
    parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
    parser.add_argument("--show_figs", nargs="?", help="True/False")
    parser.add_argument("--save_figs", nargs="?", help="True/False")
    parser.add_argument("--tag", nargs="?", help="raw/final")
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig

    if args.show_figs == 'True':
        show_figs = True
    elif args.show_figs == 'False':
        show_figs = False
    else:
        show_figs = False
        print('not able to interperet show_figs flag, must be "True" or "False", default = False')
    print('show_figs: ' + str(show_figs))

    if args.save_figs == 'True':
        save_figs = True
    elif args.save_figs == 'False':
        save_figs = False
    else:
        save_figs = False
        print('not able to interperet save_figs flag, must be "True" or "False", default = False')
    print('save_figs: ' + str(save_figs))

    if args.tag == 'raw' or 'final':
        tag = args.tag
    else:
        NameError('not able to interperet --tag, must be "raw" or "final", datatype = string')

    # hardcoded file names
    if tag == 'raw':
        experiment_file_name = 'fly.hdf5'
    else:
        experiment_file_name = 'fly_' + tag + '.hdf5'

    json_file_name = 'fly.json'
    response_set_name_prefix = 'mask_' #once channel is added, will be of form mask_ch1 (these are names saved from process_data.py)

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
    print('func_channels_num = ' + str(func_channels_num))

    #TODO: other data needed from fly_json?

    # file_path (complete path to .hdf5 file)
    experiment_file_path = os.path.join(experiment_file_directory, experiment_file_name)

    print('experiment_file_path: ' + repr(os.path.join(experiment_file_directory, experiment_file_name)))

    # instatiate correct data object plugin based on rig

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
    
        
    
    # get series numbers (series_num) from hdf5
    
    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path))) # datatype = list of strings but individual series numbers will be converted to int before using methods

    # load imaging object with first series to get fly metadata
    ID = imaging_data.ImagingDataObject(experiment_file_path,
                                        int(series_num[0]),
                                        quiet=False)
    

    # get fly metadata for later
    fly_metadata = ID.getSubjectMetadata()
    print('fly_metadata: ' + repr(fly_metadata))

    # IMPORTANT: SET TIMING_CHANNEL_IND (visanalysis assumes 0 and will default to first photodiode (ie fly left) of not set!)
    prep = fly_metadata['prep']
    
    if prep == 'fly right optic lobe': #TODO: add condition to check if there are multiple PD channels?  Currently assuming there are two PD recordings
        timing_channel_ind = 1
    elif prep == 'fly left optic lobe':
        timing_channel_ind = 0
    else:
        'could not find photodiode channel based on prep, defaulting to 0'
        timing_channel_ind = 0
    
    print('photodiode timing_channel_ind: ' + repr(timing_channel_ind))

    # initialize data structures (dicts) to store data for all series and channels
    roi_data = {}
    volume_frame_offsets = {}
    run_parameters = {}
    epoch_parameters = {}
    acquisition_metadata = {}
    unique_intensity_values = {}
    unique_center_index_values = {}
    unique_radius_values = {}
    unique_parameter_values = {}
    mean_response = {}
    sem_response = {}
    trial_response_by_stimulus = {}

    ## Extract all series and channel data into dicts for analysis

    for current_series in series_num: #loop through all series
        current_series = int(current_series) # methods expect series number to be datatype int
        sn = 'sn' + str(current_series)
        
        #update imaging object with current series number

        plug.updateImagingDataObject(experiment_file_directory, 
                                    experiment_file_name,
                                    current_series)
        
        ID = plug.ImagingDataObject
        ID.timing_channel_ind = timing_channel_ind # IMPORTANT: set timing channel index for photodiode
        ## PARAMETERS & METADATA

        # volume_frame_offsets: dict (key = sn)
        volume_frame_offsets[sn] = ID.getVolumeFrameOffsets()
        #print('volume_frame_offset keys: ' + repr(volume_frame_offsets.keys()))

        # run_parameters: dict (key = sn) of dicts (key = parameter name)
        run_parameters[sn] = ID.getRunParameters()
        #print('run_parameters keys: ' + repr(run_parameters.keys()))
        #print('run_parameters[sn] keys: ' + repr(run_parameters[sn].keys()))

        # epoch_parameters:  dict (key = sn) of dicts (key = parameter name) of all epoch parameters, one for each epoch (trial)
        epoch_parameters[sn] = ID.getEpochParameters()[0] #method returns list with one element which is dict (key=parameter names)
        #print('epoch_parameters keys: ' + repr(epoch_parameters.keys()))
        #print('epoch_parameters[sn] keys: ' + repr(epoch_parameters[sn].keys()))

        # acquisition_metadata: dict (key = sn) of dicts (key = parameter name)
        acquisition_metadata[sn] = ID.getAcquisitionMetadata()
        #print('acquisition_metadata keys: ' + repr(acquisition_metadata.keys()))
        #print('acquisition_metadata[sn] keys: ' + repr(acquisition_metadata[sn].keys()))


        for current_channel in func_channels_num: #loop through channels
            current_channel = int(current_channel) # methods expect channel number to be datatype int            
            ch = 'ch' + str(current_channel)
            response_set_name = response_set_name_prefix + ch
            
            ## Save ROIS AND RESPONSES for each series and channel to roi_data
            
            # You can access the aligned region response data just as with hand-drawn rois, using the 'aligned' prefix argument
            
            #roi_data: dict (key = sn,ch)
            #of dicts (key = roi_response/roi_image/roi_mask/epoch_response/time_vector)
            #NOTE:roi_response is a list where each entry is a roi
            roi_data[sn,ch] = ID.getRoiResponses(response_set_name, roi_prefix='aligned')

            # check stimulus timing from photodiode trace
            

            #print('roi_data keys: ' + repr(roi_data.keys()))
            #print('roi_data[sn,ch] keys: ' + repr(roi_data[sn,ch].keys()))
            print('current series, channel: ' + str(sn) + ', ' + str(ch))

            # extract unique values of all parameters

            # extract data by intensity, center location, radius
            unique_parameter_values[sn], mean_response[sn,ch], sem_response[sn,ch], trial_response_by_stimulus[sn,ch] = ID.getTrialAverages(roi_data[sn,ch]['epoch_response'], parameter_key=['angle'])
        
            print('unique_parameter_values: ' + repr(unique_parameter_values))

            # unique_parameter_values: dict (key = sn) of lists of each unique value of parameter_key (?)
            # mean_response: dict (key = sn,ch) of numpy arrays (nroi x unique values of parameter_key x time)
            # sem_response: dict (key = sn,ch) of numpy arrays (nroi x unique values of parameter_key x time)
            # trial_response_by_stimulus: dict (key = sn,ch) of lists for each unique value of parameter_key (?)

            # print('unique_parameter_values keys: ' + repr(unique_parameter_values.keys()))
            # print('unique_parameter_values datatype: ' + repr(type(unique_parameter_values[sn])))
            # print('unique_parameter_values size: ' + repr(len(unique_parameter_values[sn])))
            # print('mean_response keys: ' + repr(mean_response.keys()))
            # print('mean_response datatype: ' + repr(type(mean_response[sn,ch])))
            # print('mean_response size: ' + repr(mean_response[sn,ch].shape))
            # print('sem_response keys: ' + repr(sem_response.keys()))
            # print('sem_response datatype: ' + repr(type(sem_response[sn,ch])))
            # print('trial_response_by_stimulus keys: ' + repr(trial_response_by_stimulus.keys()))
            # print('trial_response_by_stimulus datatype: ' + repr(type(trial_response_by_stimulus[sn,ch])))
            # print('trial_response_by_stimulus size: ' + repr(len(trial_response_by_stimulus[sn,ch])))
            
    #extract number of raw rois from hdf5 (assumes all series have same number of rois which is true if process_data.py was used to make rois)
    sn = 'sn' + series_num[0] #use  first series
    ch = 'ch' + struct_channel_num[0] # use structural channel (arbitrary but shouldnt matter)
    n_roi = len(roi_data[sn,ch]['roi_response'])
    print('n_roi = ' + str(n_roi))

    ## convert mask format to bool
    roi_mask_bool = bruker.convertMaskToBool(roi_data[sn,ch]['roi_mask'])
    print('roi_mask_bool shape: ' + repr(roi_mask_bool.shape))


    """ 
    plots to make:

        selecting rois:
            1. mean responses to search stimulus - all rois on one plot
                2 rows (per condition (light/dark flash))
                2 columns (per channel)
            2. mean responses to search stimulus - individual rois
                2 rows (per condition (light/dark flash))
                2 columns (per channel)
        responses:    
            3. avg roi intensity (over imaging session), one fig for each series
            4. individual responses (short flashes)
                two channel responses for 2 series, plus avg roi intensity over entire time, 
                2 rows (per light condition (light/dark))
                2 columns (per channel)
            5. individual responses (long flashes)
                two channel responses for 2 series, plus avg roi intensity over entire time, 
                2 rows (per light condition (light/dark))
                2 columns (per channel)
            6. mean responses (short flashes)
                two channel responses for 2 series, plus avg roi intensity over entire time, 
                2 rows (per light condition (light/dark))
                2 columns (per channel)
            7. mean responses (long flashes)
                two channel responses for 2 series, plus avg roi intensity over entire time, 
                2 rows (per light condition (light/dark))
                2 columns (per channel)
        
        roi_image: 
            8. all rois overlayed on each series functional scan (structural channel)
                4 panel subplot (one per slice if 3d)
            9. each roi overlayed on each series functional scan (structural channel)
                4 panel subplot

        save: raw_rois folder
            files:
                1. search_mean_response_roi_#_
                2. mean_intensity_over_series_2_roi_#_
                3. mean_intensity_over_series_3_roi_#_
                4. individual_responses_flash_25ms_roi_#_
                5. individual_responses_flash_300ms_roi_#_
                6. mean_response_flash_25ms_roi_#_
                7. mean_response_flash_300ms_roi_#_
                8. all_rois_image_slice_#_
                9. roi_image_slice_#_


        repeat for final_rois after roi selection
    
    
    """

    # figure out which series is which

    series = 'sn' + series_num[0] #assume there is only one series
    print('analyzing series: ' + series)


    # load correct center locations if necessary

    if tag == 'final':
        print('todo: loading preferred direction for each roi')
        #TODO: load preferred direction info from json file

        # roi_centers_json_file_name = 'final_roi_centers.json'
        
        # #extract center locations for each roi
        # with open(pathlib.Path(experiment_file_directory, roi_centers_json_file_name), 'r') as file:
        #     roi_centers_final_json = json.load(file)
        # roi_centers = [int(x.split(',')[-1]) for x in roi_centers_final_json['roi_centers']]
        #  # keep second elements of list of tuples (roi_ind, center_ind)
        # print('roi_centers: ' + repr(roi_centers))

    # make raw fig save directory (if it doesn't exist)
    figs_folder_name = tag + '_roi_figs'
    figs_dir = os.path.join(experiment_file_directory,figs_folder_name)
    os.makedirs(figs_dir, exist_ok=True)

    # Plot: overlay all roi_masks on roi_image
    fig_name_string = 'mask_overlay_all'
    fig_format = '.pdf'
    ch = 'ch' + struct_channel_num[0]
    for current_series in series_num:
        sn = 'sn' + current_series
        if len(roi_mask_bool.shape)==3: #data is 2d
            im = plot_tools.overlayImage(roi_data[sn,ch]['roi_image'], roi_mask_bool, alpha=0.5, colors=None)
            fig = plt.subplot(1,1,1)
            fig.imshow(im)
            if save_figs:
                plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)
                fig_name = fig_name_string + '_{}_rois_series_{}_'.format(tag,current_series)
            if show_figs:
                plt.show()
        elif len(roi_mask_bool.shape)==4:
            for slice in range(roi_mask_bool.shape[3]):
                im = plot_tools.overlayImage(roi_data[sn,ch]['roi_image'][:,:,slice], roi_mask_bool[:,:,:,slice], alpha=0.5, colors=None)
                fig = plt.subplot(1,1,1)
                fig.imshow(im)
                if save_figs:
                    fig_name = fig_name_string + '_{}_rois_series_{}_slice_{}_'.format(tag,current_series,slice)
                    plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)
                if show_figs:
                    plt.show()

    # # Plot: individual roi_mask on roi_image (relevant slices)
    # fig_name_string = 'mask_overlay'
    # fig_format = '.pdf'
    # ch = 'ch' + struct_channel_num[0]
    # for current_series in series_num:
    #     sn = 'sn' + current_series
    #     for roi_ind in range(n_roi):
    #         if len(roi_mask_bool.shape)==3: #data is 2d
    #             im = plot_tools.overlayImage(roi_data[sn,ch]['roi_image'], roi_mask_bool[roi_ind,...], alpha=0.5, colors=None)
    #             fig = plt.subplot(1,1,1)
    #             fig.imshow(im)
    #             if save_figs:
    #                 fig_name = fig_name_string + '_roi_{}_series_{}_'.format(roi_ind,current_series)
    #                 plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)
    #             if show_figs:
    #                 plt.show()
    #         elif len(roi_mask_bool.shape)==4:
    #             for slice in range(roi_mask_bool.shape[3]):
    #                 if np.any(roi_mask_bool[roi_ind,:,:,slice]):
    #                     im = plot_tools.overlayImage(roi_data[sn,ch]['roi_image'][:,:,slice], roi_mask_bool[:,:,:,slice], alpha=0.5, colors=None)
    #                     fig = plt.subplot(1,1,1)
    #                     fig.imshow(im)
    #                     if save_figs:
    #                         fig_name = fig_name_string + '_roi_{}_series_{}_slice_{}_'.format(roi_ind,current_series,slice)
    #                         plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)
    #                     if show_figs:
    #                         plt.show()


    # plot 1: mean responses to search stimulus (series 1) - all rois plotted together
                #2 rows (per condition (light/dark flash))
                #2 columns (per channel)

    #[plot_tools.cleanAxes(x) for x in ax.ravel()]
    sn = series
    fig_format = '.pdf'
    #[x.set_ylim([-0.2, 0.2]) for x in ax.ravel()] # better way to set ax limits???  could find max of mean_responses for example
    for ch_ind, current_channel in enumerate(func_channels_num): #loop through channels   
        fig_name = 'mean_response_ch{}_all_{}_rois'.format(current_channel,tag)     
        fh, ax = plt.subplots(2, len(unique_parameter_values[sn])/2, figsize=(12, 12*(9/16)),constrained_layout = True)
        ch = 'ch' + current_channel
        for u_ind, up in enumerate(unique_parameter_values[sn]):
            ax[ch_ind, u_ind].plot(roi_data[sn,ch]['time_vector'], mean_response[sn,ch][:, u_ind, :].T)
            ax[ch_ind, u_ind].set_title('Ch{}, Angle = {}'.format(current_channel,up))
            ax[ch_ind, u_ind].set_ylabel('Mean Response (dF/F)')
            ax[ch_ind, u_ind].set_xlabel('Time (s)')
            ax[ch_ind, u_ind].axvspan(run_parameters[sn]['pre_time'], run_parameters[sn]['pre_time'] + run_parameters[sn]['stim_time'], color='gray', alpha=0.2)
    plt.suptitle("mean response, all {} rois".format(tag))

    if save_figs:
        plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)

    if show_figs:
        plt.show()

    plt.close()

    # plot 2: mean responses to stimulus (series 1) - per roi
                #2 rows (per condition (light/dark flash))
                #2 columns (per channel)
    sn = series #assume first series is search stim
    fig_name_string = 'mean_response'
    fig_format = '.pdf'
    fh, ax = plt.subplots(len(func_channels_num), len(unique_parameter_values[sn]), figsize=(12, 12*(9/16)),constrained_layout = True)
    #[plot_tools.cleanAxes(x) for x in ax.ravel()]
    for roi_ind in range(n_roi):
        #[x.set_ylim([-0.2, 0.2]) for x in ax.ravel()] # better way to set ax limits???  could find max of mean_responses for example
        for ch_ind, current_channel in enumerate(func_channels_num): #loop through channels        
            ch = 'ch' + current_channel
            for u_ind, up in enumerate(unique_intensity_values[sn]):
                ax[ch_ind, u_ind].plot(roi_data[sn,ch]['time_vector'], mean_response[sn,ch][roi_ind, u_ind, :].T)
                ax[ch_ind, u_ind].set_title('Ch{}, Angle = {}'.format(current_channel,up))
                ax[ch_ind, u_ind].set_ylabel('Mean Response (dF/F)')
                ax[ch_ind, u_ind].set_xlabel('Time (s)')
                ax[ch_ind, u_ind].axvspan(run_parameters[sn]['pre_time'], run_parameters[sn]['pre_time'] + run_parameters[sn]['stim_time'], color='gray', alpha=0.2)
        plt.suptitle("mean response, {} roi {} ".format(tag,roi_ind))
        
        if save_figs:
            fig_name = fig_name_string + '_roi_{}_'.format(roi_ind)
            plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)

        if show_figs:
            plt.show()

        plt.close()
    
    

    # plot 3: avg roi intensity over series (1, 2, 3) - for each roi
    fig_name_string = 'mean_intensity'
    fig_format = '.pdf'
    for roi_ind in range(n_roi):
        fh, ax = plt.subplots(len(func_channels_num), len(series_num), figsize=(12, 12*(9/16)),constrained_layout = True)
        for series_ind, current_series in enumerate(series_num):
            sn = 'sn' + current_series
            for ch_ind, current_channel in enumerate(func_channels_num):
                ch = 'ch' + current_channel
                ax[ch_ind, series_ind].plot(roi_data[sn,ch]['roi_response'][roi_ind].T)
                ax[ch_ind, series_ind].set_xlabel('Frame')
                ax[ch_ind, series_ind].set_ylabel('Avg ROI intensity')
                ax[ch_ind, series_ind].set_title('series {}, channel {}'.format(current_series, current_channel))
        plt.suptitle('mean roi intensity, {} roi {}'.format(tag,roi_ind))

        if save_figs:
            fig_name = fig_name_string + '_roi_{}_'.format(roi_ind)
            plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)

        if show_figs:
            plt.show()

        plt.close()

    # plot 4: avg roi intensity over series (1, 2, 3) - for all rois
    fig_name = 'mean_intensity_all_{}_rois'.format(tag)
    fig_format = '.pdf'
    fh, ax = plt.subplots(len(func_channels_num), len(series_num), figsize=(12, 12*(9/16)),constrained_layout = True)
    for series_ind, current_series in enumerate(series_num):
        sn = 'sn' + current_series
        for ch_ind, current_channel in enumerate(func_channels_num):
            ch = 'ch' + current_channel
            ax[ch_ind, series_ind].plot(np.mean(roi_data[sn,ch]['roi_response'][:], axis=0).T)
            ax[ch_ind, series_ind].set_xlabel('Frame')
            ax[ch_ind, series_ind].set_ylabel('Avg ROI intensity')
            ax[ch_ind, series_ind].set_title('series {}, channel {}'.format(current_series, current_channel))
    plt.suptitle('mean roi intensity, all {} rois'.format(tag))

    if save_figs:
        plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)

    if show_figs:
        plt.show()

    plt.close()
    

    if tag == 'good' or tag == 'final':

        # FIGURE: individual roi data for smallest radius circle at each center location (to pick correct center location)

        sn = series
        fig_format = '.pdf'
        fig_name_string = 'mean_responses_by_angle'

        for roi_ind in range(n_roi):

            fh, ax = plt.subplots(len(func_channels_num), len(unique_intensity_values[sn]), figsize=(12, 12*(9/16)),constrained_layout = True)
            # plot response for all center locations on same axes for each intensity, channel

            min_radius = min(unique_parameter_values[sn], key=lambda x: x[2])[2] # smallest radius

            for ch_ind, current_channel in enumerate(func_channels_num):
                ch = 'ch' + current_channel
                for intensity_ind, intensity in enumerate(unique_intensity_values[sn]):
                    for u_ind, up in enumerate(unique_parameter_values[sn]):
                        current_intensity = up[0]
                        current_center_index = up[1]
                        current_radius = up[2]
                        if current_intensity == intensity:
                            if current_radius == min_radius:
                                ax[ch_ind, intensity_ind].plot(roi_data[sn,ch]['time_vector'], mean_response[sn,ch][roi_ind, u_ind, :].T, label='center index: {}'.format(current_center_index))
                                ax[ch_ind, intensity_ind].legend(loc='upper right')
                                ax[ch_ind, intensity_ind].set_ylabel('Response (dF/F)')
                                ax[ch_ind, intensity_ind].set_xlabel('Time (s)')
                                ax[ch_ind, intensity_ind].set_title('Ch{}, Intensity = {}'.format(current_channel,current_intensity))

                                ax[ch_ind, intensity_ind].axvspan(run_parameters[sn]['pre_time'], run_parameters[sn]['pre_time'] + run_parameters[sn]['stim_time'], color='gray', alpha=0.2)
            plt.suptitle('Mean responses, {} roi {}, radius {}'.format(tag,roi_ind,min_radius))

            if save_figs:
                fig_name = fig_name_string + '_roi_{}_'.format(roi_ind)
                plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)

            if show_figs:
                plt.show()

            plt.close()
        
    if tag == 'final':

        # FIGURE: analyze all data from correct center location per roi

        sn = mapping_series
        fig_format = '.pdf'
        fig_name_string = 'mean_responses_by_radii'

        for roi_ind in range(n_roi):

            fh, ax = plt.subplots(len(func_channels_num), len(unique_intensity_values[sn]), figsize=(12, 12*(9/16)),constrained_layout = True)
            # plot response for all radii on same axes at correct center location for each intensity, channel

            center_index = roi_centers[roi_ind]

            for ch_ind, current_channel in enumerate(func_channels_num):
                ch = 'ch' + current_channel
                for intensity_ind, intensity in enumerate(unique_intensity_values[sn]):
                    for u_ind, up in enumerate(unique_parameter_values[sn]):
                        current_intensity = up[0]
                        current_center_index = up[1]
                        current_radius = up[2]
                        if current_intensity == intensity:
                            if current_center_index == center_index:
                                ax[ch_ind, intensity_ind].plot(roi_data[sn,ch]['time_vector'], mean_response[sn,ch][roi_ind, u_ind, :].T, label='radius: {}'.format(current_radius))
                                ax[ch_ind, intensity_ind].legend(loc='upper right')
                                ax[ch_ind, intensity_ind].set_ylabel('Response (dF/F)')
                                ax[ch_ind, intensity_ind].set_xlabel('Time (s)')
                                ax[ch_ind, intensity_ind].set_title('Ch{}, Intensity = {}'.format(current_channel,current_intensity))

                                ax[ch_ind, intensity_ind].axvspan(run_parameters[sn]['pre_time'], run_parameters[sn]['pre_time'] + run_parameters[sn]['stim_time'], color='gray', alpha=0.2)
            plt.suptitle('Mean responses, {} roi {}'.format(tag,roi_ind))

            if save_figs:
                fig_name = fig_name_string + '_roi_{}_'.format(roi_ind)
                plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)

            if show_figs:
                plt.show()

            plt.close()

        # FIGURE: analyze mean data from correct center location

        sn = mapping_series
        fig_format = '.pdf'
        fig_name_string = 'on-center_mean_responses_by_radii'

        fh, ax = plt.subplots(len(func_channels_num), len(unique_intensity_values[sn]), figsize=(12, 12*(9/16)),constrained_layout = True)
        # plot response for all radii on same axes at correct center location for each intensity, channel

        # extract on-center data for all rois

        on_center_mean_response = np.empty([n_roi, len(func_channels_num), len(unique_intensity_values[sn]), len(unique_radius_values[sn]), mean_response[sn,'ch' + str(func_channels_num[0])].shape[2]]) # numpy arrays(roi x channel x intensity x radius x time)
        on_center_sem_response = np.empty(on_center_mean_response.shape) # numpy array (roi x channel x intensity x radius x time)

        for roi_ind in range(n_roi):

            center_index = roi_centers[roi_ind]

            for ch_ind, current_channel in enumerate(func_channels_num):
                ch = 'ch' + current_channel
                for u_ind, up in enumerate(unique_parameter_values[sn]):
                    current_intensity = up[0]
                    intensity_ind = unique_intensity_values[sn].index(current_intensity)
                    current_center_index = up[1]
                    current_radius = up[2]
                    radius_ind = unique_radius_values[sn].index(current_radius)
                    if current_center_index == center_index:
                        on_center_mean_response[roi_ind, ch_ind, intensity_ind, radius_ind,:] = mean_response[sn,ch][roi_ind, u_ind, :]
                        on_center_sem_response[roi_ind, ch_ind, intensity_ind, radius_ind,:] = sem_response[sn,ch][roi_ind, u_ind, :]

        print('on_center_mean_response: ' + repr(on_center_mean_response.shape))

        # plot on-center mean responses

        for ch_ind, current_channel in enumerate(func_channels_num):
            ch = 'ch' + current_channel
            for intensity_ind, current_intensity in enumerate(unique_intensity_values[sn]):
                for radius_ind, current_radius in enumerate(unique_radius_values[sn]):
                    ax[ch_ind, intensity_ind].plot(roi_data[sn,ch]['time_vector'], np.mean(on_center_mean_response[:, ch_ind, intensity_ind, radius_ind, :], axis=0).T, label='radius: {}'.format(current_radius))
                    ax[ch_ind, intensity_ind].legend(loc='upper right')
                    ax[ch_ind, intensity_ind].set_ylabel('Response (dF/F)')
                    ax[ch_ind, intensity_ind].set_xlabel('Time (s)')
                    ax[ch_ind, intensity_ind].set_title('Ch{}, Intensity = {}'.format(current_channel,current_intensity))
            plt.suptitle('On-center mean responses')

        if save_figs:
            fig_name = fig_name_string
            plt.savefig(os.path.join(figs_dir,fig_name + fig_format), dpi=400, transparent=True)

        if show_figs:
            plt.show()

        plt.close()
