"""
attach metadata

this step comes before drawing rois and extracting data

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""

import sys
import os

from visanalysis.plugin import base as base_plugin

# call structure: sbatch visanalysis_preprocess.py --dataset_path "path" --experiment_name "name"
     #--series_number number --rig "rigID"

# dataset_path: string that contains the path to the folder that has the .hdf5 file
# experiment: experiment name
# series: series number
# rig: rigID (Bruker or AODscope)

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

