"""
analyze data:
1. open hdf5 file

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

    
    ##TODO: select rois to keep

    for current_series in series_num: #loop through all series
        current_series = int(current_series) # methods expect series number to be datatype int
        #update imaging object with current series number

        plug.updateImagingDataObject(experiment_file_directory, #NOTE: does this set current_series in imaging data object? I think so (instantiates self.ImagingObject with current_series)
                                    experiment_file_name,
                                    current_series)
        

        for current_channel in func_channels_num: #loop through channels
            current_channel = int(current_channel) # methods expect channel number to be datatype int            
            
            ID = plug.ImagingDataObject

            ID.getVolumeFrameOffsets()

            ##v PARAMETERS & METADATA

            # all run_parameters as a dict
            run_parameters = ID.getRunParameters()

            # specified run parameter
            protocol_ID = ID.getRunParameters('protocol_ID')
            print(protocol_ID)

            # epoch_parameters: list of dicts of all epoch parameters, one for each epoch (trial)
            epoch_parameters = ID.getEpochParameters()
            # Pass a param key to return a list of specified param values, one for each trial
            current_intensity = ID.getEpochParameters('current_intensity')

            # fly_metadata: dict
            fly_metadata = ID.getFlyMetadata()
            prep = ID.getFlyMetadata('prep')
            print(prep)

            # acquisition_metadata: dict
            acquisition_metadata = ID.getAcquisitionMetadata()
            sample_period = ID.getAcquisitionMetadata('sample_period')
            print(sample_period)
            
            ## ROIS AND RESPONSES

            # Get list of rois present in the hdf5 file for this series
            roi_set_names = ID.getRoiSetNames()
            print(roi_set_names)
            
            # getRoiResponses() wants a ROI set name, returns roi_data (dict)
            channel = 'ch' + str(current_channel)
            roi_set_name = 'mask_' + channel
            
            roi_data = {}
            roi_data[channel] = ID.getRoiResponses(roi_set_name)

            roi_data[channel].keys()

        # See the ROI overlaid on top of the image
        # ID.generateRoiMap(roi_name='set1', z=1)

        ## Plot whole ROI response across entire series
        fh0, ax0 = plt.subplots(1, 1, figsize=(12, 4))
        ax0.plot(roi_data.get('roi_response')[0].T)
        ax0.set_xlabel('Frame')
        ax0.set_ylabel('Avg ROI intensity')

        ## Plot ROI response for all trials
        # 'epoch_response' is shape (rois, trials, time)
        fh1, ax1 = plt.subplots(1, 1, figsize=(6, 4))
        ax1.plot(roi_data.get('time_vector'), roi_data.get('epoch_response')[0, :, :].T)
        ax1.set_ylabel('Response (dF/F)')
        ax1.set_xlabel('Time (s)')

        ## Plot trial-average responses by specified parameter name

        unique_parameter_values, mean_response, sem_response, trial_response_by_stimulus = ID.getTrialAverages(roi_data.get('epoch_response'), parameter_key='current_intensity')
        roi_data.keys()

        fh, ax = plt.subplots(1, len(unique_parameter_values), figsize=(10, 2))
        [x.set_ylim([-0.2, 0.2]) for x in ax.ravel()]
        #[plot_tools.cleanAxes(x) for x in ax.ravel()]
        for u_ind, up in enumerate(unique_parameter_values):
            ax[u_ind].plot(roi_data['time_vector'], mean_response[:, u_ind, :].T)
            ax[u_ind].set_title('Flash = {}'.format(up))
            ax[u_ind].set_ylabel('Mean Response (dF/F)')
            ax[u_ind].set_xlabel('Time (s)')
        #plot_tools.addErrorBars(ax[0], roi_data.get('time_vector'), roi_data.get('epoch_response')[0, :, :].T, stat='sem')
        # max suggests looking a function fillbetween for sem error bars
        #plot_tools.addScaleBars(ax[0], dT=2.0, dF=0.10, T_value=-0.5, F_value=-0.14)


        ## Loop through ROIs and choose which to keep

        roi_data_1 = ID.getRoiResponses('test_cells_ch1')
        unique_parameter_values_1, mean_response_1, sem_response_1, trial_response_by_stimulus_1 = ID.getTrialAverages(roi_data_1.get('epoch_response'), parameter_key='current_intensity')

        roi_data_2 = ID.getRoiResponses('test_cells_ch2')
        unique_parameter_values_2, mean_response_2, sem_response_2, trial_response_by_stimulus_2 = ID.getTrialAverages(roi_data_2.get('epoch_response'), parameter_key='current_intensity')

        unique_parameter_values = unique_parameter_values_1; # assuming same parameter values for both channels for now

        #concatenate two channels
        mean_response = np.append(np.expand_dims(mean_response_1,3), np.expand_dims(mean_response_2,3),3)
        sem_response = np.append(np.expand_dims(sem_response_1,3), np.expand_dims(sem_response_2,3),3)
        trial_response_by_stimulus = np.append(np.expand_dims(trial_response_by_stimulus_1,4), np.expand_dims(trial_response_by_stimulus_2,4),4)

        # keep_roi_ind = []

        for roi_ind in range(np.shape(mean_response)[0]):
            fh, ax = plt.subplots(2, len(unique_parameter_values), figsize=(10, 4))
            for param_ind, up in enumerate(unique_parameter_values):
                for ch_ind in [0,1]:
                    ax[param_ind, ch_ind].plot(roi_data['time_vector'], mean_response[roi_ind, param_ind, :,ch_ind])
                    ax[param_ind, ch_ind].set_title('Flash = {}, Roi = {}'.format(up,roi_ind))
                    ax[param_ind, ch_ind].set_ylabel('Mean Response (dF/F)')
                    ax[param_ind, ch_ind].set_xlabel('Time (s)')
            # fh.show()
            # response = input('Keep ROI?  y or n')
            # if response == 'y':
            #     keep_roi_ind.append(roi_ind)
            

        keep_roi_ind = [] # enter good rois manually for now

            ##TODO: resave hdf5 with selected rois only