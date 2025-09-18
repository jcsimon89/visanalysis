"""
Example script to interact with ImagingDataObject, and extract data.

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
# %%
from visanalysis.analysis import imaging_data, shared_analysis
from visanalysis.util import plot_tools
import matplotlib.pyplot as plt
import os
import pickle
import numpy as np


# %% initialize

data_directory = {} #dictionary with keys as condition names and values as paths to data directories
filenames = {} #dictionary with keys as condition names and values as lists of pkl filenames
voxel_mean = {} #dictionary with keys as condition names and values as dictionaries of channel:voxel_mean arrays
f0_data = {} #dictionary with keys as (condition name, filename) tuples and values as loaded f0_data dictionaries
driver_fly = 'JS256'
keys = ['dual sensor', 'pde', 'pde*']
effector_fly = ['JS252', 'JS257', 'JS258']

for ind, key in enumerate(keys):
    data_directory[key] = os.path.join('C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/{}_x_{}'.format(driver_fly,effector_fly[ind]))
    print('data_directory[{}]: {}'.format(key,data_directory[key]))

# %% Filter and load datafiles



for key in keys:
    filenames[key] = shared_analysis.filterPklFiles(data_directory=data_directory[key],
                    file_search_string='/fly*/f0_data.pkl',
                    quiet=False,
                    recursive=True)
    print('series[{}]: '.format(key) + str(filenames[key]))

# %% Load datafiles

for key in keys:
    for filename in filenames[key]:
        print('Loading: {}'.format(filename))
        with open(filename, 'rb') as f:
            f0_data[key,filename] = pickle.load(f)
            voxel_mean_temp = f0_data[key,filename]['voxel_mean']
        if key not in voxel_mean:
            voxel_mean[key] = voxel_mean_temp
        else:
            for ch in voxel_mean_temp.keys():
                voxel_mean[key][ch] = np.append(voxel_mean[key][ch], voxel_mean_temp[ch], axis=0)
            
        print('loaded voxel_mean from file: {}'.format(filename))
    for ch in voxel_mean[key]:
        print('voxel_mean[{}][{}] shape: {}'.format(key,ch,voxel_mean[key][ch].shape))



# %% Plot group results by condition


shared_analysis.plotF0ByConditionComparison(voxel_mean,quiet=False)

