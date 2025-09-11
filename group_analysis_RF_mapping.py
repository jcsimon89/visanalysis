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

data_directory = os.path.join('C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS256_x_JS252')


# %% Filter datafiles


mapping_series = shared_analysis.filterDataFiles(data_directory=data_directory,
                    file_search_string='/fly*/fly_final.hdf5',
                    target_fly_metadata={},
                    target_series_metadata={'stim_time': 0.3, 'radius': [5,10,15,20]},
                    exclude_series_numbers=[],
                    target_roi_series=[],
                    target_groups=[],
                    quiet=True,
                    recursive=True)



print('RF Mapping Series: ')
for series in mapping_series:
    print('file_name - ' + str(series['file_name']))
    print('series - ' + str(series['series']))

# %% Plot group results by condition

mapping_IDs=[]
for series in mapping_series:
    mapping_IDs.append(imaging_data.ImagingDataObject(series['file_name'], series['series'], quiet=True))

shared_analysis.plotAllResponsesByCondition_RF_mapping(mapping_IDs, ch_names=['mask_ch1','mask_ch2'], condition='intensity', roi_prefix='aligned')
# %% ImagingDataObject wants a path to an hdf5 file and a series number from that file
# ID = imaging_data.ImagingDataObject(file_path,
#                                     series_number,
#                                     quiet=False)

# %% Quickly look at roi responses (averaged across all trials)

# shared_analysis.plotRoiResponses(ID, roi_name='glom')

# %%

# shared_analysis.plotResponseByCondition(ID, roi_name='set_2', condition='current_intensity', eg_ind=0)