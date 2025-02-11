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
    func_channels_num = [func_channels[i].split('_')[-1] for i in range(len(func_channels))] #datatype=list of strings
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

    print('experiment_file_name: ' + repr(os.path.join(experiment_file_directory, experiment_file_name)))

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
    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path))) # datatype = list of strings
    print('series_num = '+ str(series_num))

    ## start response extraction and analysis

    for current_series in series_num: #loop through all series
        
        #update imaging object with current series number

        plug.updateImagingDataObject(experiment_file_directory,
                                    experiment_file_name,
                                    current_series)
        

        for current_channel in func_channels_num: #loop through channels (TODO: string or int? TODO: get channel_num from fly.json)
            
            #derive image file name and path

            if current_channel == '2':
                image_file_name = 'channel_2_moco_bg_func.nii'
            else:
                image_file_name = 'channel_' + current_channel + '_moco_func.nii'


            image_relative_directory = 'func' + str(int(current_series)-1) + '/moco' #folder where .nii is, assumes func_ folder counting starts from 0 which series counter starts from 1
            image_file_directory = os.path.join(experiment_file_directory, image_relative_directory)
            image_file_path = os.path.join(image_file_directory, image_file_name)
            
            #associate image data

            plug.updateImageSeries(data_directory=image_file_directory, #NOTE: doesnt do  anything with current_series!
                                    image_file_name=image_file_name,
                                    series_number=current_series,
                                    channel=current_channel)
        
            
            # Save region responses and mask to data file #TODO: make sure not to overwrite mask!

            response_set_name = 'mask_ch' + current_channel

            plug.saveRegionResponsesFromMask(file_path=experiment_file_path,
                                                series_number=series_number,
                                                response_set_name=response_set_name,
                                                mask=roi_mask,
                                                include_zero=False)

            # ID = imaging_data.ImagingDataObject(experiment_file_path,
            #                             current_series,
            #                             quiet=False)
            
            ID = plug.ImagingDataObject

            # Mask-aligned roi data gets saved under /aligned
            # Hand-drawn roi data gets saved under /rois
            print(ID.getRoiSetNames(roi_prefix='rois'))
            print(ID.getRoiSetNames(roi_prefix='aligned'))
            # You can access the aligned region response data just as with hand-drawn rois, using the 'aligned' prefix argument
            roi_data = ID.getRoiResponses(response_set_name, roi_prefix='aligned')

            #TODO: Plot region responses and masks

            #TODO: select rois to keep (based on results from which series or all?)
            #TODO: resave hdf5 with selected rois only
            #


            

# how to select rois: for each roi, could plot both channel responses for series 1, 2, 3 and then pick
# could instead just look at search stimulus