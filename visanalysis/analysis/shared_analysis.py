"""
Shared analysis tools.

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
from visanalysis.util import plot_tools
from collections.abc import Sequence
from scipy.interpolate import interp1d
import scipy.stats


def matchQuery(epoch_parameters, query):
    """

    params:
        epoch_parameters: single epoch_parameter dict
        query: dict, key-value pairs indicate matching parameter editable_values
            e.g. query = {'current_intensity': 0.75}

    Returns:
        Bool, True if all param values match, false otherwise
    """
    return np.all([epoch_parameters.get(key) == query[key] for key in query])


def filterTrials(epoch_response, ID, query, return_inds=False):
    no_trials = epoch_response.shape[1]
    matching_trials = np.where([matchQuery(ep, query) for ep in ID.getEpochParameters()[:no_trials]])[0]

    if return_inds:
        return epoch_response[:, matching_trials, :], matching_trials
    else:
        return epoch_response[:, matching_trials, :]


def getUniqueParameterCombinations(param_keys, ID):
    ep_params = [[ep.get(x, None) for x in param_keys]for ep in ID.getEpochParameters()]
    return list({tuple(row) for row in ep_params})


def plotAllResponsesByCondition(ImagingDataObjects, ch_names, condition, roi_prefix='rois'):
    # plot all roi responses by 
    # unpack roi_data and unique_parameter_values from all ImagingDataObjects

    #ImagingDataObjects is a list of ImagingDataObject instances
    #ch_names is a list of roi names (ex ['mask_ch1', 'mask_ch2'])

    roi_data={}
    unique_parameter_values=[]
    mean_response={}
    sem_response={}
    trial_response_by_stimulus={}
    sample_period = []
    n_timepoints = []
    n_roi_total = 0
    for exp_ind, ImagingData in enumerate(ImagingDataObjects):
        for ch_ind, ch_name in enumerate(ch_names):
            # get roi data
            roi_data[exp_ind,ch_ind] = ImagingData.getRoiResponses(ch_name, roi_prefix=roi_prefix)
            # extract mean_response and unique_parameter_values by condition 
            unique_parameter_values, mean_response[exp_ind,ch_ind], sem_response[exp_ind,ch_ind], trial_response_by_stimulus[exp_ind,ch_ind] = ImagingData.getTrialAverages(roi_data[exp_ind,ch_ind]['epoch_response'], parameter_key='intensity')

            #print('roi_data["epoch_response"].shape: {}'.format(roi_data['epoch_response'].shape))
            #print('roi_data["roi_response"][0].shape: {}'.format(roi_data['roi_response'][0].shape))
            #print('roi_data["time_vector"].shape: {}'.format(roi_data['time_vector'].shape))
            #roi_data['epoch_response'] - 3D array of responses for each trial (roi, trial, time)
            #roi_data['roi_response'] - list of roi responses (roi, time)
            #roi_data['time_vector'] - 1d array of timepoints (one set for all measurements)
            #print('type(roi_data): {}'.format(type(roi_data)))
            #print('roi_data.keys(): {}'.format(roi_data.keys()))
            #print('.time_vector: {}'.format(roi_data.get('time_vector')))
            n_roi = mean_response[exp_ind,ch_ind].shape[0]
            n_roi_total = n_roi_total + n_roi
            n_timepoints.append(len(roi_data[exp_ind,ch_ind]['time_vector']))
            sample_period.append(ImagingData.getAcquisitionMetadata('sample_period'))
        run_parameters = ImagingData.getRunParameters() # will be redefined in loop, but should be same for all scans

    # need to enforce same sample period for all experiments, assume for now
  
    # combine data intelligently into single array with single time vector
    # interpolate to match longest time vector
    print('unique_parameter_values: ' + repr(unique_parameter_values))
    n_timepoints = max(n_timepoints)
    sample_period = min(sample_period)
    total_time = n_timepoints * sample_period
    print('resampling to {} timepoints, {}s sample period, {}s per epoch'.format(n_timepoints, sample_period, total_time))
    # assume same unique_parameter_values for all experiments
    print('unique_parameter_values ({}): {}'.format(condition,unique_parameter_values))

    frames= range(0,n_timepoints)
    time_vector = frames*sample_period
    # mean_responses[exp_ind,ch_ind].shape: (nroi x unique values of parameter_key x time)

    #mean_response(nroi x unique values of parameter_key x time)
    #mean_response_interp(nroi x unique values of parameter_key x time x ch)
    #mean_response_interp = np.empty([n_roi_total,len(unique_parameter_values),n_timepoints,len(ch_names)])
    #print('mean_response_interp.shape: ' + repr(mean_response_interp.shape))

    for exp_ind, ImagingData in enumerate(ImagingDataObjects):
        for ch_ind, ch_name in enumerate(ch_names):
            # interpolate
            #f = interp1d(roi_data[exp_ind,roi_ind]['time_vector'],roi_data[exp_ind,roi_ind]['epoch_response'],kind='linear',axis=2)
            f = interp1d(roi_data[exp_ind,ch_ind]['time_vector'],mean_response[exp_ind,ch_ind][:,:,:],kind='linear',axis=2,bounds_error = False)
            if ch_ind==0:
                response_interp_temp = np.expand_dims(f(time_vector),axis=-1) # add channel dim
            else:
                response_interp_temp = np.append(response_interp_temp, np.expand_dims(f(time_vector),axis=-1),axis=-1) # add channel dim, expand along ch axis
            print('response_interp_temp.shape: ' + repr(response_interp_temp.shape))
            #mean_response_interp[:,:,:,ch_ind] = response_interp_temp
        if exp_ind==0:
                # first assignment
                mean_response_interp = response_interp_temp
        else:
                mean_response_interp = np.append(mean_response_interp,response_interp_temp, axis=0)

    print('mean_response_interp.shape: ' + repr(mean_response_interp.shape))
    fh1, ax1 = plt.subplots(len(ch_names), len(unique_parameter_values), figsize=(10, 10*9/16),constrained_layout = True)
    fh2, ax2 = plt.subplots(len(ch_names), len(unique_parameter_values), figsize=(10, 10*9/16),constrained_layout = True)
    for u_ind, u_value in enumerate(unique_parameter_values):
        for ch_ind, ch_name in enumerate(ch_names):
            if 'ch1' in ch_name:
                ch_label = 'ch1'
            elif 'ch2' in ch_name:
                ch_label = 'ch2'
            else:
                ch_label = 'unk ch'
                print('could not extract channel label form roi set name')
            #query = {condition: u_value}
            #trials = filterTrials(roi_data.get('epoch_response'), ImagingData, query)
            y = np.mean(mean_response_interp[:,u_ind,:,ch_ind],axis=0).T
            error = scipy.stats.sem(mean_response_interp[:,u_ind,:,ch_ind],axis=0).T
            ax1[ch_ind,u_ind].plot(time_vector, mean_response_interp[:,u_ind,:,ch_ind].T)
            ax1[ch_ind,u_ind].set_title('{}, Intensity = {}, {}ms Flash'.format(ch_label,u_value,1000*run_parameters['stim_time'])) #, linestyle='-', color=ImagingData.colors[0])
            ax1[ch_ind, u_ind].set_ylabel('Response (dF/F)')
            ax1[ch_ind, u_ind].set_xlabel('Time (s)')
            ax1[ch_ind, u_ind].axvspan(run_parameters['pre_time'], run_parameters['pre_time'] + run_parameters['stim_time'], color='gray', alpha=0.2)
            ax2[ch_ind,u_ind].plot(time_vector,y, color = 'k')
            ax2[ch_ind,u_ind].fill_between(time_vector, y-error, y+error, color = 'c', alpha = 0.3) #, linestyle='-', color=ImagingData.colors[0])
            ax2[ch_ind,u_ind].set_title('{}, Intensity = {}, {}ms Flash'.format(ch_label,u_value,1000*run_parameters['stim_time']))
            ax2[ch_ind, u_ind].set_ylabel('Response (dF/F)')
            ax2[ch_ind, u_ind].set_xlabel('Time (s)')
            ax2[ch_ind, u_ind].axvspan(run_parameters['pre_time'], run_parameters['pre_time'] + run_parameters['stim_time'], color='gray', alpha=0.2)

def plotRoiResponses(ImagingData, roi_name):
    roi_data = ImagingData.getRoiResponses(roi_name)

    fh, ax = plt.subplots(1, int(roi_data.get('epoch_response').shape[0]+1), figsize=(6, 2))
    [x.set_axis_off() for x in ax]
    [x.set_ylim([-0.25, 1]) for x in ax]

    for r_ind in range(roi_data.get('epoch_response').shape[0]):
        time_vector = roi_data.get('time_vector')
        no_trials = roi_data.get('epoch_response')[r_ind, :, :].shape[0]
        current_mean = np.mean(roi_data.get('epoch_response')[r_ind, :, :], axis=0)
        current_std = np.std(roi_data.get('epoch_response')[r_ind, :, :], axis=0)
        current_sem = current_std / np.sqrt(no_trials)

        ax[r_ind].plot(time_vector, current_mean, 'k')
        ax[r_ind].fill_between(time_vector,
                               current_mean - current_sem,
                               current_mean + current_sem,
                               alpha=0.5)
        ax[r_ind].set_title(int(r_ind))

        if r_ind == 0:  # scale bar
            plot_tools.addScaleBars(ax[r_ind], 1, 1, F_value=-0.1, T_value=-0.2)


def filterDataFiles(data_directory,
                    file_search_string='*.hdf5',
                    target_fly_metadata={},
                    target_series_metadata={},
                    exclude_series_numbers=[],
                    target_roi_series=[],
                    target_groups=[],
                    quiet=False,
                    recursive=False):
    """
    Searches through a directory of visprotocol datafiles and finds datafiles/series that match the search values
    Can search based on any number of fly metadata params or run parameters

    Params
        -data_directory: directory of visprotocol data files to search through
        -target_fly_metadata: (dict) key-value pairs of target parameters to search for in the fly metadata
        -target_series_metadata: (dict) key-value pairs of target parameters to search for in the series run (run parameters)
        -target_roi_series: (list) required roi_series names
        -target_groups: (list) required names of groups under series group

    Returns
        -matching_series: List of matching series dicts with all fly & run params as well as file name and series number
    """
    fileNames = glob.glob(data_directory + file_search_string, recursive=recursive)
    if not quiet:
        print('filename search string: {}'.format(file_search_string))
        print('searching in directory: {}'.format(data_directory))
        print('recursive search: {}'.format(recursive))
        print('Found {} files in {}'.format(len(fileNames), data_directory))
        print('fileNames: ' + repr(fileNames))

    # collect key/value pairs for all series in data directory
    all_series = []
    for ind, fn in enumerate(fileNames):

        with h5py.File(fn, 'r') as data_file:
            for fly in data_file.get('Subjects'):
                fly_metadata = {}
                for f_key in data_file.get('Subjects').get(fly).attrs.keys():
                    fly_metadata[f_key] = data_file.get('Subjects').get(fly).attrs[f_key]

                for epoch_run in data_file.get('Subjects').get(fly).get('epoch_runs'):
                    series_metadata = {}
                    for s_key in data_file.get('Subjects').get(fly).get('epoch_runs').get(epoch_run).attrs.keys():
                        series_metadata[s_key] = data_file.get('Subjects').get(fly).get('epoch_runs').get(epoch_run).attrs[s_key]

                    new_series = {**fly_metadata, **series_metadata}
                    new_series['series'] = int(epoch_run.split('_')[1])
                    new_series['file_name'] = fn

                    existing_roi_sets = list(data_file.get('Subjects').get(fly).get('epoch_runs').get(epoch_run).get('rois').keys())
                    new_series['rois'] = existing_roi_sets
                    existing_groups = list(data_file.get('Subjects').get(fly).get('epoch_runs').get(epoch_run).keys())
                    new_series['groups'] = existing_groups

                    all_series.append(new_series)

    # search in all series for target key/value pairs
    match_dict = {**target_fly_metadata, **target_series_metadata}
    matching_series = []
    for series in all_series:
        if checkAgainstTargetDict(match_dict, series):
            if np.all([r in series.get('rois') for r in target_roi_series]):
                if np.all([r in series.get('groups') for r in target_groups]):
                    matching_series.append(series)

    matching_series = [series for series in matching_series if series.get('series') not in exclude_series_numbers]

    matching_series = sorted(matching_series, key=lambda d: d['file_name'] + '-' + str(d['series']).zfill(3))

    # filter by series number

    if not quiet:
        print('Found {} matching series'.format(len(matching_series)))
    return matching_series


def checkAgainstTargetDict(target_dict, test_dict):
    for key in target_dict:
        if key in test_dict:
            if not areValsTheSame(target_dict[key], test_dict[key]):
                return False  # Different values
        else:
            return False  # Target key not in this series at all

    return True


def areValsTheSame(target_val, test_val):

    if isinstance(target_val, str):
        return target_val.casefold() == test_val.casefold()
    elif isinstance(target_val, bool):
        if isinstance(test_val, str):
            return str(target_val).casefold() == test_val.casefold()

        return target_val == test_val

    elif isinstance(target_val, (int, float)):  # Scalar
        if isinstance(test_val, (int, float)):
            return float(target_val) == float(test_val)  # Ignore type for int vs. float here
        else:
            return False
    elif isinstance(target_val, (Sequence, np.ndarray)):  # Note already excluded possibility of string by if ordering
        if isinstance(test_val, (Sequence, np.ndarray)):
            # Ignore order of arrays, and ignore float vs. int
            return np.all(np.sort(target_val, axis=0).astype(float) == np.sort(test_val, axis=0).astype(float))
        else:
            return False

    else:
        print('----')
        print('Unable to match ')
        print('Target {} ({})'.format(target_val, type(target_val)))
        print('Test {} ({})'.format(test_val, type(test_val)))
        print('----')
        return False
