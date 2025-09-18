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

def plotF0ByConditionComparison(search_IDs):
    """
    Plot average F0 by condition for multiple groups of ImagingDataObjects.

    Parameters
    ----------
    search_IDs : dict
        Dictionary with keys as group names and values as lists of ImagingDataObjects.
    ch_names : list, optional
        List of channel names to plot. The default is ['mask_ch1','mask_ch2'].
    condition : str, optional
        Condition to group by. The default is 'intensity'.
    roi_prefix : str, optional
        Prefix for ROI names to select. The default is 'aligned'.

    Returns
    -------
    None.

    """
    
    group_data = {}
    
    for key, id_list in search_IDs.items():
        group_data[key] = []
        for ID in id_list:
            roi_names = [roi for roi in ID.getRoiNames() if roi.startswith(roi_prefix)]
            for roi in roi_names:
                data = ID.getRoiResponses(roi_name=roi, ch_names=ch_names, response_type='dF/F0')
                metadata = ID.getSeriesMetadata()
                for cond_value in set(metadata[condition]):
                    cond_indices = [i for i, val in enumerate(metadata[condition]) if val == cond_value]
                    avg_response = np.mean([data[i] for i in cond_indices], axis=0)
                    group_data[key].append((cond_value, avg_response))
    
    plt.figure(figsize=(10, 6))
    
    for key, responses in group_data.items():
        cond_values = sorted(set([resp[0] for resp in responses]))
        avg_responses = [np.mean([resp[1] for resp in responses if resp[0] == cv], axis=0) for cv in cond_values]
        
        for ch_idx, ch_name in enumerate(ch_names):
            plt.subplot(len(ch_names), 1, ch_idx + 1)
            for cv_idx, cv in enumerate(cond_values):
                plt.plot(avg_responses[cv_idx], label=f'{key} - {condition}: {cv}')
            plt.title(f'Channel: {ch_name}')
            plt.xlabel('Time')
            plt.ylabel('dF/F0')
            plt.legend()
    
    plt.tight_layout()
    plt.show()

# %% initialize

data_directory = {}
series = {}
series_IDs = {}
driver_fly = 'JS140'
keys = ['dual sensor', 'pde', 'pde*']
effector_fly = ['JS252', 'JS257', 'JS258']

for ind, key in enumerate(keys):
    data_directory[key] = os.path.join('C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/{}_x_{}'.format(driver_fly,effector_fly[ind]))
    print('data_directory[{}]: {}'.format(key,data_directory[key]))

# %% Filter datafiles

for key in data_directory.keys():
    series[key] = shared_analysis.filterDataFiles(data_directory=data_directory[key],
                    file_search_string='/fly*/f0_data.hdf5',
                    target_fly_metadata={},
                    target_series_metadata={},
                    exclude_series_numbers=[],
                    target_roi_series=[],
                    target_groups=[],
                    quiet=True,
                    recursive=True)

    print('series[{}]: '.format(key) + str(series[key]))

for key in series.keys():
    current_series = series[key]
    print('{} series: '.format(key))
    print('file_name - ' + str(current_series['file_name']))
    print('series - ' + str(current_series['series']))

    for series in series[key]:
        if key not in series_IDs:
            series_IDs[key] = []
        series_IDs[key].append(imaging_data.ImagingDataObject(series['file_name'], series['series'], quiet=True))

# %% Plot group results by condition





plotF0ByConditionComparison(series_IDs)

