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
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from visanalysis.plugin import base as base_plugin

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

    structural_channel = fly_json['structural_channel']
    structural_channel_num = structural_channel.split('_')[-1]

    #TODO: other data needed from fly_json?

    # derive image_file_path
    # assumes folder structure from snake_eyesss

        # assuming structural channel is brightest functional channel and will be used for roi selection,
        # currently only doing bg subtraction for channel 2

    if structural_channel_num == '2':
        image_file_name = 'channel_2_moco_bg_func.nii'
    else:
        image_file_name = 'channel_' + structural_channel_num + '_moco_func.nii'

    image_directory = 'func' + str(int(series_number)-1) + '/moco' #folder where .nii is, assumes func_ folder counting starts from 0 which series counter starts from 1
    image_file_path = os.path.join(experiment_file_directory, image_directory, image_file_name)
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

    for series in range(number_of_series): #loop through all series
    
        #update imaging object with current series number

        plug.updateImagingDataObject(experiment_file_path,
                                    experiment_file_name,
                                    series)

        for channel in channels: #loop through channels
        
            #associate raw image data
            
            plug.updateImageSeries(experiment_file_path,
                                    image_file_name='TSeries-20210707-001_reg.nii',
                                    series_number=series_number,
                                    channel=0)
        
            
            
            
            # Save region responses and mask to data file #TODO: make sure not to overwrite mask!
            plug.saveRegionResponsesFromMask(file_path=file_path,
                                                series_number=series_number,
                                                response_set_name='mask_1',
                                                mask=mask,
                                                include_zero=False)

            # Retrieve saved mask region responses from data file

            ID = imaging_data.ImagingDataObject(file_path,
                                        series_number,
                                        quiet=False)

            # Mask-aligned roi data gets saved under /aligned
            # Hand-drawn roi data gets saved under /rois
            ID.getRoiSetNames(roi_prefix='aligned')

            # You can access the aligned region response data just as with hand-drawn rois, using the 'aligned' prefix argument
            roi_data = ID.getRoiResponses('mask_1', roi_prefix='aligned')

            # Plot region responses and masks
            z_slice = 2
            fh, ax = plt.subplots(2, 2, figsize=(8, 3),
                                gridspec_kw={'width_ratios': [1, 4]})
            [x.set_axis_off() for x in ax.ravel()]

            colors = 'rgb'
            for r_ind, response in enumerate(roi_data['roi_response']):
                ax[r_ind, 0].imshow((roi_data['roi_mask'][:, :, z_slice] == (r_ind+1)).T)
                ax[r_ind, 1].plot(response, color=colors[r_ind])

