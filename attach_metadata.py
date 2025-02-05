"""
attach metadata

this step comes before drawing rois and extracting data

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""

import sys
import os
import argparse

from visanalysis.plugin import base as base_plugin

# call structure: python visanalysis_preprocess.py --dataset_path "path" --experiment_name "name" --rig "rigID"

# dataset_path: string that contains the path to the folder that has the .hdf5 file
# experiment: experiment name
# rig: rigID (Bruker or AODscope)


if __name__ == '__main__':
    # parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset_path", nargs="?", help="Folder pointing to hdf5")
    parser.add_argument("--experiment_name", nargs="?", help="hdf5 file name")
    parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
    args = parser.parse_args()

    dataset_path = args.dataset_path
    experiment_name = args.experiment_name
    rig = args.rig

    experiment_file_path = os.path.join(dataset_path, experiment_name + '.hdf5')

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

    plug.attachData(experiment_name, dataset_path)

    print('Attached data to {}'.format(experiment_name))

