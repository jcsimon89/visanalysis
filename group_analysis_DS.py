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

data_directory = os.path.join('C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS144_x_JS252')

# %% Filter datafiles


series = shared_analysis.filterDataFiles(data_directory=data_directory,
                    file_search_string='/fly*/fly_final.hdf5',
                    target_fly_metadata={},
                    target_series_metadata={'stim_time': 4},
                    exclude_series_numbers=[],
                    target_roi_series=[],
                    target_groups=[],
                    quiet=True,
                    recursive=True)



print('DS Series: ')
for current_series in series:
    print('file_name - ' + str(current_series['file_name']))
    print('series - ' + str(current_series['series']))

# %% Plot group results by condition

IDs=[]
for current_series in series:
    IDs.append(imaging_data.ImagingDataObject(current_series['file_name'], current_series['series'], quiet=True))

shared_analysis.plotAllResponses(IDs, ch_names=['Lobula_ch1','Lobula_ch2'], roi_prefix='aligned')
# %% ImagingDataObject wants a path to an hdf5 file and a series number from that file
# ID = imaging_data.ImagingDataObject(file_path,
#                                     series_number,
#                                     quiet=False)

# %% Quickly look at roi responses (averaged across all trials)

# shared_analysis.plotRoiResponses(ID, roi_name='glom')

# %%

# shared_analysis.plotResponseByCondition(ID, roi_name='set_2', condition='current_intensity', eg_ind=0)