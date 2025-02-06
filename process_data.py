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

# call structure: python process_data.py --data_directory "path" --rig "rigID"

# data_directory: string that contains the path to the folder that has the .hdf5 file
# rig: rigID (Bruker or AODscope)

#NOTE: all file names and paths are relative to the data_directory which is the main fly folder ie /fly_001


if __name__ == '__main__':
    ## parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_directory", nargs="?", help="Folder pointing to hdf5")
    parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
    args = parser.parse_args()

    data_directory = args.data_directory
    rig = args.rig

    # hardcoded file names
    experiment_file_name = 'fly.hdf5'
    json_file_name = 'fly.json'
    
    #extract info from fly.json 

    with open(pathlib.Path(data_directory, json_file_name), 'r') as file:
        fly_json = json.load(file)

    structural_channel = fly_json['structural_channel']
    structural_channel_num = structural_channel.split('_')[-1]

    #TODO: other data needed from fly_json?

    #TODO extract subject number from hdf5 or folder name? not sure if needed

    print('experiment_file_name: ' + repr(os.path.join(data_directory, experiment_file_name)))

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

    plug.attachData(experiment_file_name, data_directory)

    print('Attached metadata to {}'.format(experiment_file_name))

    ##draw roi masks using reduced GUI
    
    gui_path = str(os.path.join(os.getcwd(),"gui/DataGUI_prog.py"))

    os.system('python ' + gui_path
              + ' --data_directory ' + data_directory
              + ' --experiment_file_name ' + experiment_file_name
              + ' --rig ' + rig)

    # series_number = 1

    # plug.updateImagingDataObject(data_directory,
    #                             experiment_file_name,
    #                             series_number)

    
    
    # plug.updateImageSeries(data_directory,
    #                         image_file_name='TSeries-20210707-001_reg.nii',
    #                         series_number=series_number,
    #                         channel=0)
    

