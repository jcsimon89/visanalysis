"""
process data:
1. attach metadata
2. open hdf5 file

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

# call structure: python process_data.py --experiment_file_directory "path" --rig "rigID" --series_number "series_number"

# experiment_file_directory: string that contains the path to the folder that has the .hdf5 file
# rig: rigID (Bruker or AODscope)

#NOTE: all file names and paths are relative to the experiment_file_directory which is the main fly folder ie /fly_001


if __name__ == '__main__':
    ## parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment_file_directory", nargs="?", help="Folder pointing to hdf5")
    parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
    parser.add_argument("--series_number", nargs="?", help="number as string for series to use for roi selection")
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig
    series_number = args.series_number

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

    # derive image_file_path
    # assumes folder structure from snake_eyesss

        # assuming structural channel is brightest functional channel and will be used for roi selection,
        # currently only doing bg subtraction for channel 2

    if struct_channel_num[0] == '2':
        image_file_name = 'channel_2_moco_bg_func.nii'
    else:
        image_file_name = 'channel_' + struct_channel_num[0] + '_moco_func.nii'

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

    plug.attachData(experiment_file_name, experiment_file_directory)

    print('Attached metadata to {}'.format(experiment_file_name))

    ##draw roi masks using reduced GUI
    
    gui_path = str(os.path.join(os.getcwd(),"gui/DataGUI_prog.py"))

    os.system('python ' + gui_path
              + ' --experiment_file_directory ' + experiment_file_directory
              + ' --experiment_file_name ' + experiment_file_name
              + ' --experiment_file_path ' + experiment_file_path
              + ' --rig ' + rig
              + ' --series_number ' + series_number
              + ' --image_file_path ' + image_file_path)

    
    
    # Retrieve saved mask region responses from data file

    ID = imaging_data.ImagingDataObject(experiment_file_path,
                                        series_number,
                                        quiet=False)

    roi_data = ID.getRoiMasks('roi_set_name')
    roi_mask_binary = roi_data['roi_mask'] #shape: roi_index, x, y ,(z)
    roi_image = roi_data['roi_image'] #shape: x,y,(z)

    print('dimensions of roi_image: ' + str(np.shape(roi_image)))
    print('dimensions of roi_mask_binary: ' + str(np.shape(roi_mask_binary)))

    # figure out if data is volume or slice
    data_spacial_dim = len(roi_image.shape)
   
    # convert roi mask to numbers for each separate roi
    roi_mask = np.zeros(roi_mask_binary.shape[1:])
    print('dimensions of roi_mask: ' + str(np.shape(roi_mask))) # shape:x,y,(z)

    if data_spacial_dim == 3: #data is a volume
        for roi_ind in range(roi_mask_binary.shape[0]):
            for i in range(roi_mask_binary.shape[1]):
                for j in range(roi_mask_binary.shape[2]):
                    for k in range(roi_mask_binary.shape[3]):
                        if roi_mask_binary[roi_ind,i,j,k] == True:
                            roi_mask[i,j,k]=roi_ind+1 #since roi_ind starts at 0 and we want to label the first roi with value=1

    elif data_spacial_dim == 2: #data is a slice
        for roi_ind in range(roi_mask_binary.shape[0]):
            for i in range(roi_mask_binary.shape[1]):
                for j in range(roi_mask_binary.shape[2]):
                    if roi_mask_binary[roi_ind,i,j] == True:
                            roi_mask[i,j]=roi_ind+1 #since roi_ind starts at 0 and we want to label the first roi with value=1
    else:
        print('data does not appear to be a volume or a slice according to the dimensions of roi_image')

    print('unique values in roi_mask: ' + str(np.unique(roi_mask)))

    ## extract other important metadata for analysis
    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path))) # datatype = list of strings but individual series numbers will be converted to int before using methods
    print('series_num = '+ str(series_num))

    ## start response extraction and analysis

    for current_series in series_num: #loop through all series
        current_series = int(current_series) # methods expect series number to be datatype int
        #update imaging object with current series number

        plug.updateImagingDataObject(experiment_file_directory, #NOTE: does this set current_series in imaging data object? I think so (instantiates self.ImagingObject with current_series)
                                    experiment_file_name,
                                    current_series)
        

        for current_channel in func_channels_num: #loop through channels (TODO: string or int? TODO: get channel_num from fly.json)
            current_channel = int(current_channel) # methods expect channel number to be datatype int
            #derive image file name and path

            if current_channel == 2:  # assuming we only do background subtraction on channel 2 (green)
                image_file_name = 'channel_' + str(current_channel) + '_moco_bg_func.nii' 
            else:
                image_file_name = 'channel_' + str(current_channel) + '_moco_func.nii' 


            image_relative_directory = 'func' + str(int(current_series)-1) + '/moco' #folder where .nii is, assumes func_ folder counting starts from 0 which series counter starts from 1
            image_file_directory = os.path.join(experiment_file_directory, image_relative_directory)
            image_file_path = os.path.join(image_file_directory, image_file_name)
            
            #associate image data

            plug.updateImageSeries(data_directory=image_file_directory, #NOTE: doesnt do  anything with current_series!
                                    image_file_name=image_file_name,
                                    series_number=current_series,
                                    channel=current_channel)
        
            
            # Save region responses and mask to data file #TODO: make sure not to overwrite mask!

            response_set_name = 'mask_ch' + str(current_channel)

            plug.saveRegionResponsesFromMask(file_path=experiment_file_path,
                                                series_number=current_series,
                                                response_set_name=response_set_name,
                                                mask=roi_mask,
                                                include_zero=False)

            ID = imaging_data.ImagingDataObject(experiment_file_path,
                                        current_series,
                                        quiet=False)

            # Mask-aligned roi data gets saved under /aligned
            # Hand-drawn roi data gets saved under /rois

            ID.getRoiSetNames(roi_prefix='aligned')

            # You can access the aligned region response data just as with hand-drawn rois, using the 'aligned' prefix argument
            roi_data = ID.getRoiResponses(response_set_name, roi_prefix='aligned')

            #TODO: Plot region responses and masks, save fig
           
           
            z_slice = 2
            fh, ax = plt.subplots(len(roi_data), 2, figsize=(8, 3),
                      gridspec_kw={'width_ratios': [1, 4]})
            [x.set_axis_off() for x in ax.ravel()]

            colors = 'rgb'
            for r_ind, response in enumerate(roi_data['roi_response']):
                print('r_ind = ' + str(r_ind))
                ax[r_ind, 0].imshow((roi_mask[:, :, z_slice] == (r_ind+1)).T)
                ax[r_ind, 1].plot(response, color=colors[r_ind])
            
##TODO: select rois to keep (based on results from which series or all?)

for current_series in series_num: #loop through all series
        current_series = int(current_series) # methods expect series number to be datatype int
        #update imaging object with current series number

        plug.updateImagingDataObject(experiment_file_directory, #NOTE: does this set current_series in imaging data object? I think so (instantiates self.ImagingObject with current_series)
                                    experiment_file_name,
                                    current_series)
        

        for current_channel in func_channels_num: #loop through channels (TODO: string or int? TODO: get channel_num from fly.json)
            current_channel = int(current_channel) # methods expect channel number to be datatype int            
            
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
            roi_data = ID.getRoiResponses('2_9')#,background_subtraction=True) # option for background subtraction, requires roi named 'bg'
            #roi_data = ID.getRoiResponses('AllMedullaGreen')
            roi_data.keys()

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