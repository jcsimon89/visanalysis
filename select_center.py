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
    parser.add_argument("--input_tag", nargs="?", help="Good/Final/etc. label for input data (can be empty string)")
    parser.add_argument("--output_tag", nargs="?", help="Good/Final/etc. label for output data")
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig
    
    if args.input_tag is None:
        input_tag = ''
    else:
        input_tag = args.input_tag
        
    output_tag = args.output_tag
    print('input_tag:' + input_tag)

    if args.save == 'True':
        save = True
    elif args.save == 'False':
        save = False
    else:
        save = False
        print('not able to interperet save flag, must be "True" or "False", default = False')
    print('save: ' + str(save))

    # file names
    json_file_name = 'fly.json'
    
    if input_tag == '':
        experiment_file_name = 'fly.hdf5'
    else:
        experiment_file_name = 'fly_' + input_tag + '.hdf5'

    roi_centers_json_file_name = output_tag + '_roi_centers.json'
    response_set_name_prefix = 'mask_' #once channel is added, will be of form mask_ch1 (these are names saved from process_data.py)
    roi_prefix = 'aligned'
    hdf5_save_path = os.path.join(experiment_file_directory,'fly_' + output_tag + '.hdf5') # save path
    if save:
        print('hdf5_save_path: ' + hdf5_save_path)

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
    
    #extract center locations for each roi from final_roi_centers.json

    with open(pathlib.Path(experiment_file_directory, roi_centers_json_file_name), 'r') as file:
        roi_centers_final_json = json.load(file)
    roi_centers = [int(x.split(',')[-1]) for x in roi_centers_final_json['roi_centers']]
    roi_ind_final = [int(x.split(',')[0]) for x in roi_centers_final_json['roi_centers']]
        # keep second elements of list of tuples (roi_ind, center_ind)
    print('roi_centers: ' + repr(roi_centers))

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

            #figure out mask
            mask = roi_data[sn,ch]['roi_mask'] #nd array with same dims as brain containing diff integer for each roi
            for i in range(mask.shape[0]):
                for j in range(mask.shape[1]):
                    try:
                        for k in range(mask.shape[2]):
                            if mask[i,j,k] != 0:
                                if mask[i,j,k] not in roi_ind_final:
                                    mask[i,j,k] = 0
                                else: mask[i,j,k] = roi_ind_final.index(mask[i,j,k]) + 1 #find ind in roi_ind_final and set to that value +1 (since mask counting starts at 1 not 0)
                                #TODO: start mask value back at 1 for final rois or keep same roi numbers?
                                # currently renumbering final rois
                    except:
                        print('data is not three dimensional, trying with two dimensions')
                        if mask[i,j] != 0:
                                if mask[i,j] not in roi_ind_final:
                                    mask[i,j] = 0
                                else: mask[i,j] = roi_ind_final.index(mask[i,j]) + 1 #find ind in roi_ind_final and set to that value +1 (since mask counting starts at 1 not 0)
                                #TODO: start mask value back at 1 for final rois or keep same roi numbers?
                                # currently renumbering final rois
            roi_data_final[sn,ch]['roi_mask'] = mask
            roi_data_final[sn,ch]['roi_image'] = roi_data[sn,ch]['roi_image'] #TODO:i think roi_image is just a meanbrain image but need to check, different for each channel???
            roi_data_final[sn,ch]['roi_response'] = [roi_data[sn,ch]['roi_response'][i] for i in roi_ind_final]

    if save:
        with h5py.File(experiment_file_path, 'r') as h5r:
            with h5py.File(hdf5_save_path, 'w') as h5w:
                
                for obj in h5r.keys():        
                    h5r.copy(obj, h5w)
                
                for roi_ind in roi_ind_final:
                    center_index = roi_centers[roi_ind_final.index(roi_ind)]

                    for epoch in h5w['/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + '/epochs']:
                                # add epoch parameter on_center (0 or 1) if center_index is roi_center

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
                        h5w.create_dataset('/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        '/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_image',data=roi_data_final[sn,ch]['roi_image'])
                        h5w.create_dataset('/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        '/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_mask',data=roi_data_final[sn,ch]['roi_mask'])
                        h5w.create_dataset('/Subjects/{}/epoch_runs/series_00'.format(subject_number) + current_series + 
                        '/{}/{}'.format(roi_prefix,response_set_name_prefix) + ch + '/roi_response',data=roi_data_final[sn,ch]['roi_response'])  
                            
        
        h5r.close()