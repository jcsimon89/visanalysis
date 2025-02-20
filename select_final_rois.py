"""
select final rois:
1. select raw_rois to keep

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

plt.ioff() #don't show plots unless called for with plt.show() 

# call structure: python analyze_data.py --experiment_file_directory "path" --rig "rigID" --save "True/False"

if __name__ == '__main__':
    ## parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment_file_directory", nargs="?", help="Folder pointing to hdf5")
    parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
    parser.add_argument("--save", nargs="?", help="True/False")
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig

    if args.save == 'True':
        save = True
    elif args.save == 'False':
        save = False
    else:
        save = False
        print('not able to interperet save flag, must be "True" or "False", default = False')
    print('save: ' + str(save))

    # hardcoded file names
    experiment_file_name = 'fly.hdf5'
    json_file_name = 'fly.json'
    final_rois_json_file_name = 'final_rois.json'
    response_set_name_prefix = 'mask_' #once channel is added, will be of form mask_ch1 (these are names saved from process_data.py)
    roi_prefix = 'aligned'
    final_hdf5_save_path = os.path.join(experiment_file_directory,'fly_final.hdf5') # save path
    if save:
        print('final_hdf5_save_path: ' + final_hdf5_save_path)

    experiment_file_path = os.path.join(experiment_file_directory,experiment_file_name)

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

    # extract roi_ind_final from final_rois.json 

    with open(pathlib.Path(experiment_file_directory, final_rois_json_file_name), 'r') as file:
        roi_ind_final_json = json.load(file)
    roi_ind_final = roi_ind_final_json['roi_ind_final'] #list of integers
    print('roi_ind_final: ' + repr(roi_ind_final))

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
   
    ## load roi_data from raw hdf5, resave data for selected rois in roi_data_final

    # load imaging object with first series to get fly metadata
    ID = imaging_data.ImagingDataObject(experiment_file_path,
                                        int(series_num[0]),
                                        quiet=False)
    
    # get fly metadata for later
    fly_metadata = ID.getSubjectMetadata()
    print('fly_metadata: ' + repr(fly_metadata))

    subject_number = fly_metadata['subject_id']

    # initialize data structures (dicts) to store data for all series and channels
    roi_data = {}

    ## Extract all series and channel data into dicts for analysis

    for current_series in series_num: #loop through all series
        current_series = int(current_series) # methods expect series number to be datatype int
        sn = 'sn' + str(current_series)
        
        #update imaging object with current series number

        plug.updateImagingDataObject(experiment_file_directory, 
                                    experiment_file_name,
                                    current_series)
        
        ID = plug.ImagingDataObject

        for current_channel in func_channels_num: #loop through channels
            current_channel = int(current_channel) # methods expect channel number to be datatype int            
            ch = 'ch' + str(current_channel)
            response_set_name = response_set_name_prefix + ch
            
            ## Save ROIS AND RESPONSES for each series and channel to roi_data
            
            # You can access the aligned region response data just as with hand-drawn rois, using the 'aligned' prefix argument
            
            #roi_data: dict (key = sn,ch)
            #of dicts (key = roi_response/roi_image/roi_mask)
            #NOTE:roi_response is a list where each entry is a roi
            roi_data[sn,ch] = ID.getRoiResponses(roi_set_name=response_set_name, roi_prefix=roi_prefix) 
            #print('roi_data keys: ' + repr(roi_data.keys()))
            #print('roi_data[sn,ch] keys: ' + repr(roi_data[sn,ch].keys()))
            #print('len(roi_data[sn,ch][roi_response]): ' + repr(len(roi_data[sn,ch]['roi_response'])))
            #print('roi_data[sn,ch][roi_image].shape: ' + repr(roi_data[sn,ch]['roi_image'].shape))
            #print('roi_data[sn,ch][roi_mask].shape: ' + repr(roi_data[sn,ch]['roi_mask'].shape))

    ##select rois to keep
    roi_data_final = {}

    for current_series in series_num:
        sn = 'sn' + current_series
        for current_channel in func_channels_num:
            ch = 'ch' + current_channel
            roi_data_final[sn,ch] = {}
            roi_data_final[sn,ch]['roi_image'] = [roi_data[sn,ch]['roi_image'][i] for i in roi_ind_final]
            roi_data_final[sn,ch]['roi_mask'] = [roi_data[sn,ch]['roi_mask'][i] for i in roi_ind_final]
            roi_data_final[sn,ch]['roi_response'] = [roi_data[sn,ch]['roi_response'][i] for i in roi_ind_final]

    if save:
        with h5py.File(experiment_file_path, 'r') as h5r:
            with h5py.File(final_hdf5_save_path, 'w') as h5w:
                for obj in h5r.keys():        
                    h5r.copy(obj, h5w)
                for current_series in series_num: #loop through all series
                    sn = 'sn' + current_series
                    for current_channel in func_channels_num: #loop through channels
                        ch = 'ch' + current_channel
                        del h5w['/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        '/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_image']
                        del h5w['/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        '/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_mask']
                        del h5w['/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        '/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_response']
                        #TODO: figure out how to save roi_mask, roi_image for final rois
                        #h5w.create_dataset('/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        #'/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_image',data=roi_data_final[sn,ch]['roi_image'])
                        #h5w.create_dataset('/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        #'/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_mask',data=roi_data_final[sn,ch]['roi_mask'])
                        h5w.create_dataset('/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        '/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_response',data=roi_data_final[sn,ch]['roi_response'])       
        
        h5r.close()