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

data_directory = os.path.join('C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS140_x_JS252')


# %% Filter datafiles

search_series = shared_analysis.filterDataFiles(data_directory=data_directory,
                    file_search_string='/fly*/fly_final.hdf5',
                    target_fly_metadata={},
                    target_series_metadata={'stim_time': 0.3, 'width_height': (25,25)},
                    exclude_series_numbers=[],
                    target_roi_series=[],
                    target_groups=[],
                    quiet=True,
                    recursive=True)

long_flash_series = shared_analysis.filterDataFiles(data_directory=data_directory,
                    file_search_string='/fly*/fly_final.hdf5',
                    target_fly_metadata={},
                    target_series_metadata={'stim_time': 0.3, 'width_height': (240,240)},
                    exclude_series_numbers=[4,5],
                    target_roi_series=[],
                    target_groups=[],
                    quiet=True,
                    recursive=True)

short_flash_series = shared_analysis.filterDataFiles(data_directory=data_directory,
                    file_search_string='/fly*/fly_final.hdf5',
                    target_fly_metadata={},
                    target_series_metadata={'stim_time': 0.025, 'width_height': (240,240)},
                    exclude_series_numbers=[],
                    target_roi_series=[],
                    target_groups=[],
                    quiet=True,
                    recursive=True)

print('Search Series: ')
for series in search_series:
    print('file_name - ' + str(series['file_name']))
    print('series - ' + str(series['series']))

print('Long Flash Series: ')
for series in long_flash_series:
    print('file_name - ' + str(series['file_name']))
    print('series - ' + str(series['series']))

print('Short Flash Series: ')
for series in short_flash_series:
    print('file_name - ' + str(series['file_name']))
    print('series - ' + str(series['series']))

# %% Plot group results by condition
search_IDs=[]
for series in search_series:
    search_IDs.append(imaging_data.ImagingDataObject(series['file_name'], series['series'], quiet=True))

shared_analysis.plotAllResponsesByCondition(search_IDs, ch_names=['mask_ch1','mask_ch2'], condition='intensity', roi_prefix='aligned')

long_flash_IDs=[]
for series in long_flash_series:
    long_flash_IDs.append(imaging_data.ImagingDataObject(series['file_name'], series['series'], quiet=True))

shared_analysis.plotAllResponsesByCondition(long_flash_IDs, ch_names=['mask_ch1','mask_ch2'], condition='intensity', roi_prefix='aligned')

short_flash_IDs=[]
for series in short_flash_series:
    short_flash_IDs.append(imaging_data.ImagingDataObject(series['file_name'], series['series'], quiet=True))

shared_analysis.plotAllResponsesByCondition(short_flash_IDs, ch_names=['mask_ch1','mask_ch2'], condition='intensity', roi_prefix='aligned')
# %% ImagingDataObject wants a path to an hdf5 file and a series number from that file
# ID = imaging_data.ImagingDataObject(file_path,
#                                     series_number,
#                                     quiet=False)

# %% Quickly look at roi responses (averaged across all trials)

# shared_analysis.plotRoiResponses(ID, roi_name='glom')

# %%

# shared_analysis.plotResponseByCondition(ID, roi_name='set_2', condition='current_intensity', eg_ind=0)