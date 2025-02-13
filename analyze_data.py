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
import h5py

# call structure: python analyze_data.py --experiment_file_directory "path" --rig "rigID"

# experiment_file_directory: string that contains the path to the folder that has the .hdf5 file

if __name__ == '__main__':
    ## parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment_file_directory", nargs="?", help="Folder pointing to hdf5")
    parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig

    # hardcoded file names
    experiment_file_name = 'fly.hdf5'
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
    print('func_channel_num = ' + str(func_channels_num))

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
    
    # save fly metadata for later
    fly_metadata = ID.getSubjectMetadata()
    print('fly_metadata: ' + repr(fly_metadata))

    # initialize data structures (dicts) to store data for all series and channels
    roi_data = {}
    volume_frame_offsets = {}
    run_parameters = {}
    epoch_parameters = {}
    acquisition_metadata = {}
    unique_intensity_values = {}
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
            #print('roi_data keys: ' + repr(roi_data.keys()))
            #print('roi_data[sn,ch] keys: ' + repr(roi_data[sn,ch].keys()))
            
            unique_intensity_values[sn], mean_response[sn,ch], sem_response[sn,ch], trial_response_by_stimulus[sn,ch] = ...
            ID.getTrialAverages(roi_data[sn,ch]['epoch_response'], parameter_key='intensity')            
            # unique_intensity_values: dict (key = sn) of dicts? (key = parameter name)
            # mean_response: dict (key = sn,ch) of dicts? (key = parameter name)
            # sem_response: dict (key = sn,ch) of dicts? (key = parameter name)
            # trial_response_by_stimulus: dict (key = sn,ch) of dicts? (key = parameter name)

            print('unique_intensity_values keys: ' + repr(unique_intensity_values.keys()))
            print('unique_intensity_values[sn] keys: ' + repr(unique_intensity_values[sn].keys()))
            print('mean_response keys: ' + repr(mean_response.keys()))
            print('mean_response[sn,ch] keys: ' + repr(mean_response[sn,ch].keys()))
            print('sem_response keys: ' + repr(sem_response.keys()))
            print('sem_response[sn,ch] keys: ' + repr(sem_response[sn,ch].keys()))
            print('trial_response_by_stimulus keys: ' + repr(trial_response_by_stimulus.keys()))
            print('trial_response_by_stimulus[sn,ch] keys: ' + repr(trial_response_by_stimulus[sn,ch].keys()))



    #TODO: Plot region responses and masks, save figs

    """ 
    plots to make:

        selecting rois:
            1. mean responses to search stimulus
                2 rows (per condition (light/dark flash))
                2 columns (per channel)
        responses:    
            2. and 3. avg roi intensity (over imaging session)
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

    # plot 1: mean responses to search stimulus (series 1)
                #2 rows (per condition (light/dark flash))
                #2 columns (per channel)
    







    ## misc plotting notes


    # See the ROI overlaid on top of the image
    # ID.generateRoiMap(roi_name='set1', z=1)

    ## Plot whole ROI response across entire series - for first roi
    fh0, ax0 = plt.subplots(1, 1, figsize=(12, 4))
    ax0.plot(roi_data[sn,ch]['roi_response'][0].T)
    ax0.set_xlabel('Frame')
    ax0.set_ylabel('Avg ROI intensity')
    plt.show()

    ## Plot ROI response for all trials
    # 'epoch_response' is shape (rois, trials, time)
    fh1, ax1 = plt.subplots(1, 1, figsize=(6, 4))
    ax1.plot(roi_data[sn,ch]['time_vector'], roi_data[sn,ch]['epoch_response'][0, :, :].T)
    ax1.set_ylabel('Response (dF/F)')
    ax1.set_xlabel('Time (s)')
    plt.show()

    ## Plot trial-average responses by specified parameter name

    unique_parameter_values, mean_response, sem_response, trial_response_by_stimulus = ID.getTrialAverages(roi_data[sn,ch]['epoch_response'], parameter_key='intensity')
    roi_data[sn,ch].keys()

    fh, ax = plt.subplots(1, len(unique_parameter_values), figsize=(10, 2))
    [x.set_ylim([-0.2, 0.2]) for x in ax.ravel()]
    #[plot_tools.cleanAxes(x) for x in ax.ravel()]
    for u_ind, up in enumerate(unique_parameter_values):
        ax[u_ind].plot(roi_data[sn,ch]['time_vector'], mean_response[:, u_ind, :].T)
        ax[u_ind].set_title('Flash = {}'.format(up))
        ax[u_ind].set_ylabel('Mean Response (dF/F)')
        ax[u_ind].set_xlabel('Time (s)')
    plt.show()
    #plot_tools.addErrorBars(ax[0], roi_data.get('time_vector'), roi_data.get('epoch_response')[0, :, :].T, stat='sem')
    # max suggests looking a function fillbetween for sem error bars
    #plot_tools.addScaleBars(ax[0], dT=2.0, dF=0.10, T_value=-0.5, F_value=-0.14)

    ## Loop through ROIs and choose which to keep

    # keep_roi_ind = []

    for roi_ind in range(np.shape(mean_response)[0]):
        fh, ax = plt.subplots(2, len(unique_parameter_values), figsize=(10, 4))
        for param_ind, up in enumerate(unique_parameter_values):
            for ch_ind in [0,1]:
                ch = 'ch' + str(ch_ind + 1)
                ax[param_ind, ch_ind].plot(roi_data['time_vector'], mean_response[roi_ind, param_ind, :,ch_ind])
                ax[param_ind, ch_ind].set_title('Flash = {}, Roi = {}'.format(up,roi_ind))
                ax[param_ind, ch_ind].set_ylabel('Mean Response (dF/F)')
                ax[param_ind, ch_ind].set_xlabel('Time (s)')
        # fh.show()
        # response = input('Keep ROI?  y or n')
        # if response == 'y':
        #     keep_roi_ind.append(roi_ind)
    

    ##TODO: select rois to keep
    keep_roi_ind = [] # enter good rois manually for now
    roi_data_final = roi_data[keep_roi_ind]

    ##TODO: more plotting, figure saving

    ##TODO: resave hdf5 with selected rois only

    h5r=h5py.File(experiment_file_path, 'r')
    final_hdf5_save_path = os.path.join(experiment_file_directory + 'fly_final.hdf5')
    
    with h5py.File(final_hdf5_save_path, 'w') as h5w:
        for obj in h5r.keys():        
            h5r.copy(obj, h5w)
        for current_series in series_num: #loop through all series
            sn = 'sn' + current_series
            for current_channel in func_channels_num: #loop through channels
                ch = 'ch' + current_channel
                h5w.create_dataset('/Subjects/1/epoch_runs/series_00' + current_series + 
                '/aligned/mask_ch' + current_channel + '/roi_response',data=roi_data_final[sn,ch]['roi_response'])       
    
    h5r.close()

#path to data in h5 file: /Subjects/1/epoch_runs/series_003/aligned/mask_ch1/roi_response
    