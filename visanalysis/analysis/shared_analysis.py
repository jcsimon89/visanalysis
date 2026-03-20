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
import seaborn as sns


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


def plotAllResponsesByCondition(ImagingDataObjects, ch_names, condition, bin_frequency, roi_prefix='rois', dff='pre'):
    """
    Plot all ROI responses by condition, pooling across ImagingDataObjects.

    Uses per-epoch time vectors and binning (instead of interpolation) to
    combine data from experiments with different sample rates or epoch lengths.

    Params:
        ImagingDataObjects: list of ImagingDataObject instances
        ch_names: list of roi set names (e.g. ['mask_ch1', 'mask_ch2'])
        condition: str, parameter key to split conditions on (e.g. 'intensity')
        bin_frequency: float, Hz. Temporal frequency for the common time grid.
        roi_prefix: str, prefix for ROI group in HDF5 file
        dff: str, dF/F method passed to getRoiResponses ('pre', 'mean', 'none')
    """
    import warnings

    response_ylabel = 'F' if dff == 'none' else r'$\Delta F/F_0$'
    bin_width = 1.0 / bin_frequency

    # ---- Step 1: Collect data from all ImagingDataObjects ----
    roi_data = {}
    unique_parameter_values = None
    epoch_inds_by_exp = {}

    for exp_ind, ImagingData in enumerate(ImagingDataObjects):
        fly_metadata = ImagingData.getSubjectMetadata()
        prep = fly_metadata.get('prep', '')
        if prep == 'fly right optic lobe':
            ImagingData.timing_channel_ind = 1
        elif prep == 'fly left optic lobe':
            ImagingData.timing_channel_ind = 0
        else:
            ImagingData.timing_channel_ind = 0
        print('Exp {}: prep={}, timing_channel_ind={}'.format(exp_ind, prep, ImagingData.timing_channel_ind))

        for ch_ind, ch_name in enumerate(ch_names):
            roi_data[exp_ind, ch_ind] = ImagingData.getRoiResponses(ch_name, roi_prefix=roi_prefix, dff=dff)

        # Get epoch groupings (same for all channels within an experiment)
        upv, ei = ImagingData.getEpochGroupingsByParameters(parameter_key=condition)
        epoch_inds_by_exp[exp_ind] = ei
        if unique_parameter_values is None:
            unique_parameter_values = upv
        run_parameters = ImagingData.getRunParameters()

    n_conditions = len(unique_parameter_values)
    n_channels = len(ch_names)

    # ---- Step 2: For each condition, pool trials across experiments ----
    # Determine global max time for bin edges
    all_max_times = []
    for exp_ind in range(len(ImagingDataObjects)):
        for tv in roi_data[exp_ind, 0]['time_vector_by_epoch']:
            if len(tv) > 0:
                all_max_times.append(tv[-1])
    global_max_time = max(all_max_times)
    bin_edges = np.arange(0, global_max_time + bin_width, bin_width)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    n_bins = len(bin_centers)

    # For each condition x channel: collect binned responses across all experiments and ROIs
    # Result shape per condition x channel: (n_rois_total, n_bins)
    binned_by_condition = {}  # key: (cond_ind, ch_ind), value: list of 1D arrays (one per roi)

    for cond_ind in range(n_conditions):
        for ch_ind in range(n_channels):
            binned_by_condition[cond_ind, ch_ind] = []

    for exp_ind in range(len(ImagingDataObjects)):
        epoch_response = roi_data[exp_ind, 0]['epoch_response']  # (rois, epochs, max_time)
        time_vector_by_epoch = roi_data[exp_ind, 0]['time_vector_by_epoch']
        n_rois = epoch_response.shape[0]

        for cond_ind in range(n_conditions):
            trial_inds = epoch_inds_by_exp[exp_ind][cond_ind]
            # Clamp to available trials
            trial_inds = trial_inds[trial_inds < epoch_response.shape[1]]
            if len(trial_inds) == 0:
                continue

            for ch_ind in range(n_channels):
                ch_epoch_response = roi_data[exp_ind, ch_ind]['epoch_response']
                ch_tvbe = roi_data[exp_ind, ch_ind]['time_vector_by_epoch']

                for roi_ind in range(ch_epoch_response.shape[0]):
                    # Bin this ROI's matching trials, then average across trials per bin
                    trial_binned = np.full((len(trial_inds), n_bins), np.nan)
                    for t_i, trial_idx in enumerate(trial_inds):
                        tv = ch_tvbe[trial_idx]
                        n_valid = len(tv)
                        if n_valid == 0:
                            continue
                        bin_indices = np.clip(np.digitize(tv, bin_edges) - 1, 0, n_bins - 1)
                        for b in range(n_bins):
                            in_bin = np.where(bin_indices == b)[0]
                            if len(in_bin) > 0:
                                trial_binned[t_i, b] = np.nanmean(ch_epoch_response[roi_ind, trial_idx, in_bin])

                    # Average across trials for this ROI -> 1D (n_bins,)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        roi_mean = np.nanmean(trial_binned, axis=0)
                    binned_by_condition[cond_ind, ch_ind].append(roi_mean)

    # Stack into arrays: (n_rois_total, n_bins)
    for cond_ind in range(n_conditions):
        for ch_ind in range(n_channels):
            if len(binned_by_condition[cond_ind, ch_ind]) > 0:
                binned_by_condition[cond_ind, ch_ind] = np.vstack(binned_by_condition[cond_ind, ch_ind])
            else:
                binned_by_condition[cond_ind, ch_ind] = np.empty((0, n_bins))

    # ---- Step 3: Plot ----
    plt.rc('font', size=16)

    # Figure 1: individual ROI traces
    fh1, ax1 = plt.subplots(n_channels, n_conditions, figsize=(10, 10 * 9 / 16), constrained_layout=True, squeeze=False)
    # Figure 2: mean ± SEM per channel
    fh2, ax2 = plt.subplots(n_channels, n_conditions, figsize=(10, 10 * 9 / 16), constrained_layout=True, squeeze=False)
    # Figure 3: mean ± SEM both channels overlaid
    fh3, ax3 = plt.subplots(1, n_conditions, figsize=(10, 5 * 9 / 16), constrained_layout=True, squeeze=False)

    ch_colors = []
    ch_labels = []
    for ch_name in ch_names:
        if 'ch1' in ch_name:
            ch_colors.append('c')
            ch_labels.append('ch1')
        elif 'ch2' in ch_name:
            ch_colors.append('orange')
            ch_labels.append('ch2')
        else:
            ch_colors.append('gray')
            ch_labels.append('unk ch')

    for cond_ind, cond_value in enumerate(unique_parameter_values):
        for ch_ind in range(n_channels):
            data = binned_by_condition[cond_ind, ch_ind]  # (n_rois_total, n_bins)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                y = np.nanmean(data, axis=0)
                error = scipy.stats.sem(data, axis=0, nan_policy='omit')

            # Fig 1: individual traces
            ax1[ch_ind, cond_ind].plot(bin_centers, data.T, alpha=0.5)
            ax1[ch_ind, cond_ind].set_title('{}, {} = {}'.format(ch_labels[ch_ind], condition, cond_value))
            ax1[ch_ind, cond_ind].set_ylabel(response_ylabel)
            ax1[ch_ind, cond_ind].set_xlabel('Time (s)')
            ax1[ch_ind, cond_ind].axvspan(run_parameters['pre_time'],
                                          run_parameters['pre_time'] + run_parameters['stim_time'],
                                          color='gray', alpha=0.2)

            # Fig 2: mean ± SEM
            ax2[ch_ind, cond_ind].plot(bin_centers, y, color='k')
            ax2[ch_ind, cond_ind].fill_between(bin_centers, y - error, y + error,
                                               color=ch_colors[ch_ind], alpha=0.4)
            ax2[ch_ind, cond_ind].set_title('{}, {} = {}'.format(ch_labels[ch_ind], condition, cond_value))
            ax2[ch_ind, cond_ind].set_ylabel(response_ylabel)
            ax2[ch_ind, cond_ind].set_xlabel('Time (s)')
            ax2[ch_ind, cond_ind].axvspan(run_parameters['pre_time'],
                                          run_parameters['pre_time'] + run_parameters['stim_time'],
                                          color='gray', alpha=0.2)

        # Fig 3: both channels overlaid
        for ch_ind in range(n_channels):
            data = binned_by_condition[cond_ind, ch_ind]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                y = np.nanmean(data, axis=0)
                error = scipy.stats.sem(data, axis=0, nan_policy='omit')
            ax3[0, cond_ind].plot(bin_centers, y, color='k')
            ax3[0, cond_ind].fill_between(bin_centers, y - error, y + error,
                                          color=ch_colors[ch_ind], alpha=0.4)
        ax3[0, cond_ind].set_title('{} = {}'.format(condition, cond_value))
        ax3[0, cond_ind].set_ylabel(response_ylabel)
        ax3[0, cond_ind].set_xlabel('Time (s)')
        ax3[0, cond_ind].axvspan(run_parameters['pre_time'],
                                 run_parameters['pre_time'] + run_parameters['stim_time'],
                                 color='gray', alpha=0.2)


def plotAllResponses(ImagingDataObjects, ch_names, roi_prefix='rois'):
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
        fly_metadata = ImagingData.getSubjectMetadata()
        print('fly_metadata: ' + repr(fly_metadata))

        # IMPORTANT: SET TIMING_CHANNEL_IND (visanalysis assumes 0 and will default to first photodiode (ie fly left) of not set!)
        prep = fly_metadata['prep']
    
        if prep == 'fly right optic lobe': #TODO: add condition to check if there are multiple PD channels?  Currently assuming there are two PD recordings
            timing_channel_ind = 1
        elif prep == 'fly left optic lobe':
            timing_channel_ind = 0
        else:
            'could not find photodiode channel based on prep, defaulting to 0'
            timing_channel_ind = 0
        print('photodiode timing_channel_ind: ' + repr(timing_channel_ind))
        ImagingData.timing_channel_ind = timing_channel_ind # IMPORTANT: set timing channel index for photodiode
        for ch_ind, ch_name in enumerate(ch_names):
            # get roi data
            roi_data[exp_ind,ch_ind] = ImagingData.getRoiResponses(ch_name, roi_prefix=roi_prefix)
            # extract mean_response and unique_parameter_values by condition 
            unique_parameter_values, mean_response[exp_ind,ch_ind], sem_response[exp_ind,ch_ind], trial_response_by_stimulus[exp_ind,ch_ind] = ImagingData.getTrialAverages(roi_data[exp_ind,ch_ind]['epoch_response'])

            #print('roi_data["epoch_response"].shape: {}'.format(roi_data['epoch_response'].shape))
            #print('roi_data["roi_response"][0].shape: {}'.format(roi_data['roi_response'][0].shape))
            #print('roi_data["time_vector"].shape: {}'.format(roi_data['time_vector'].shape))
            #roi_data['epoch_response'] - 3D array of responses for each trial (roi, trial, time)
            #roi_data['roi_response'] - list of roi responses (roi, time)
            #roi_data['time_vector'] - 1d array of timepoints (one set for all measurements)
            #print('type(roi_data): {}'.format(type(roi_data)))
            #print('roi_data.keys(): {}'.format(roi_data.keys()))
            #print('.time_vector: {}'.format(roi_data.get('time_vector')))
            if ch_ind==0: #only need to do this once for first channel
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
    #print('unique_parameter_values ({}): {}'.format(condition,unique_parameter_values))

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
    plt.rc('font', size=14)
    fh1, ax1 = plt.subplots(len(ch_names), len(unique_parameter_values), figsize=(5, 10*9/16),constrained_layout = True)
    fh2, ax2 = plt.subplots(len(ch_names), len(unique_parameter_values), figsize=(5, 10*9/16),constrained_layout = True)
    fh3, ax3 = plt.subplots(1, len(unique_parameter_values), figsize=(5, 5*9/16),constrained_layout = True)
    for ch_ind, ch_name in enumerate(ch_names):
        if 'ch1' in ch_name:
            ch_label = 'ch1'
            current_color = 'c'
        elif 'ch2' in ch_name:
            ch_label = 'ch2'
            current_color = 'orange'
        else:
            ch_label = 'unk ch'
            print('could not extract channel label form roi set name')
        #query = {condition: u_value}
        #trials = filterTrials(roi_data.get('epoch_response'), ImagingData, query)
        y = np.mean(mean_response_interp[:,0,:,ch_ind],axis=0).T
        error = scipy.stats.sem(mean_response_interp[:,0,:,ch_ind],axis=0).T
        ax1[ch_ind].plot(time_vector, mean_response_interp[:,0,:,ch_ind].T)
        ax1[ch_ind].set_title('{}, , {}ms Flash'.format(ch_label,1000*run_parameters['stim_time'])) #, linestyle='-', color=ImagingData.colors[0])
        ax1[ch_ind].set_ylabel('Response (dF/F)')
        ax1[ch_ind].set_xlabel('Time (s)')
        ax1[ch_ind].axvspan(run_parameters['pre_time'], run_parameters['pre_time'] + run_parameters['stim_time'], color='gray', alpha=0.2)
        ax2[ch_ind].plot(time_vector,y, color = 'k')
        ax2[ch_ind].fill_between(time_vector, y-error, y+error, color = current_color, alpha = 0.4) #, linestyle='-', color=ImagingData.colors[0])
        ax2[ch_ind].set_title('{}, {}ms Flash'.format(ch_label,1000*run_parameters['stim_time']))
        ax2[ch_ind].set_ylabel('Response (dF/F)')
        ax2[ch_ind].set_xlabel('Time (s)')
        ax2[ch_ind].axvspan(run_parameters['pre_time'], run_parameters['pre_time'] + run_parameters['stim_time'], color='gray', alpha=0.2)


        #query = {condition: u_value}
        #trials = filterTrials(roi_data.get('epoch_response'), ImagingData, query)
        y = np.mean(mean_response_interp[:,0,:,:],axis=0).T
        error = scipy.stats.sem(mean_response_interp[:,0,:,:],axis=0).T
        ax3.plot(time_vector,y.T, color = 'k')
        ax3.fill_between(time_vector, (y[0,:]-error[0,:]).T, (y[0,:]+error[0,:]).T, color = 'c', alpha = 0.4) #, linestyle='-', color=ImagingData.colors[0])
        ax3.fill_between(time_vector, (y[1,:]-error[1,:]).T, (y[1,:]+error[1,:]).T, color = 'orange', alpha = 0.4) #, linestyle='-', color=ImagingData.colors[0])
        ax3.set_title('{}, {}ms Flash'.format(ch_label,1000*run_parameters['stim_time']))
        ax3.set_ylabel('Response (dF/F)')
        ax3.set_xlabel('Time (s)')
        ax3.axvspan(run_parameters['pre_time'], run_parameters['pre_time'] + run_parameters['stim_time'], color='gray', alpha=0.2)

def plotAllResponsesByCondition_RF_mapping(ImagingDataObjects, ch_names, response_set_name, condition, roi_prefix='rois'):
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
    func_channels_num = len(ch_names)

    for exp_ind, ImagingData in enumerate(ImagingDataObjects):
        fly_metadata = ImagingData.getSubjectMetadata()
        print('fly_metadata: ' + repr(fly_metadata))

        # IMPORTANT: SET TIMING_CHANNEL_IND (visanalysis assumes 0 and will default to first photodiode (ie fly left) of not set!)
        prep = fly_metadata['prep']
    
        if prep == 'fly right optic lobe': #TODO: add condition to check if there are multiple PD channels?  Currently assuming there are two PD recordings
            timing_channel_ind = 1
        elif prep == 'fly left optic lobe':
            timing_channel_ind = 0
        else:
            'could not find photodiode channel based on prep, defaulting to 0'
            timing_channel_ind = 0
        print('photodiode timing_channel_ind: ' + repr(timing_channel_ind))
        ImagingData.timing_channel_ind = timing_channel_ind # IMPORTANT: set timing channel index for photodiode
        for ch_ind, ch_name in enumerate(ch_names):
            # get roi data
            roi_data[exp_ind,ch_ind] = ImagingData.getRoiResponses(ch_name, roi_prefix=roi_prefix)
            # extract mean_response and unique_parameter_values by condition 
            
            unique_intensity_values = ImagingData.getTrialAverages(roi_data[exp_ind,ch_ind]['epoch_response'], parameter_key='intensity')[0]
            unique_radius_values = ImagingData.getTrialAverages(roi_data[exp_ind,ch_ind]['epoch_response'], parameter_key='radius')[0]
            unique_center_index_values = ImagingData.getTrialAverages(roi_data[exp_ind,ch_ind]['epoch_response'], parameter_key='center_index')[0]
            unique_parameter_values, mean_response[exp_ind,ch_ind], sem_response[exp_ind,ch_ind], trial_response_by_stimulus[exp_ind,ch_ind] = ImagingData.getTrialAverages(roi_data[exp_ind,ch_ind]['epoch_response'], parameter_key=['intensity','center_index','radius'])


            #print('roi_data["epoch_response"].shape: {}'.format(roi_data['epoch_response'].shape))
            #print('roi_data["roi_response"][0].shape: {}'.format(roi_data['roi_response'][0].shape))
            #print('roi_data["time_vector"].shape: {}'.format(roi_data['time_vector'].shape))
            #roi_data['epoch_response'] - 3D array of responses for each trial (roi, trial, time)
            #roi_data['roi_response'] - list of roi responses (roi, time)
            #roi_data['time_vector'] - 1d array of timepoints (one set for all measurements)
            #print('type(roi_data): {}'.format(type(roi_data)))
            #print('roi_data.keys(): {}'.format(roi_data.keys()))
            #print('.time_vector: {}'.format(roi_data.get('time_vector')))
            if ch_ind==0: #only need to do this once for first channel
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

    on_center_mean_response_interp = np.empty((n_roi_total, func_channels_num, len(unique_intensity_values), len(unique_radius_values), n_timepoints)) # numpy arrays(roi x channel x intensity x radius x time)
    roi_counter = 0
    for exp_ind, ImagingData in enumerate(ImagingDataObjects):
        
        # aggregate on-center mean data for each roi

        roi_centers = ImagingData.getRoiParameters(roi_set_name=response_set_name, parameter='center', roi_prefix='aligned')
        print('roi_centers: ' + repr(roi_centers))


        n_roi = mean_response[exp_ind,0].shape[0] #n_roi in experiment, using first channel as reference
        for roi_ind in range(n_roi): 
            center_index = roi_centers[roi_ind]
            for ch_ind, ch_name in enumerate(ch_names):
                for u_ind, up in enumerate(unique_parameter_values):
                    current_intensity = up[0]
                    intensity_ind = unique_intensity_values.index(current_intensity)
                    current_center_index = up[1]
                    current_radius = up[2]
                    radius_ind = unique_radius_values.index(current_radius)
                    if current_center_index == center_index:
                       
                       # interpolate mean_response and sem_response to common time vector

                        f = interp1d(roi_data[exp_ind,ch_ind]['time_vector'],mean_response[exp_ind,ch_ind][roi_ind,u_ind,:],kind='linear',axis=0,bounds_error = False)

                        on_center_mean_response_interp[roi_ind + roi_counter, ch_ind, intensity_ind, radius_ind,:] = f(time_vector)
        roi_counter = roi_counter + n_roi

    print('on_center_mean_response_interp.shape (roi, channel, intensity, radius, time): ' + repr(on_center_mean_response_interp.shape))

    fig_format = '.pdf'
    fig_name_string = 'on-center_mean_responses_by_radii'

    fh1, ax1 = plt.subplots(len(ch_names), len(unique_intensity_values), figsize=(10, 10*9/16),constrained_layout = True)
    for u_ind, u_value in enumerate(unique_intensity_values):
        for ch_ind, ch_name in enumerate(ch_names):
            if 'ch1' in ch_name:
                ch_label = 'ch1'
            elif 'ch2' in ch_name:
                ch_label = 'ch2'
            else:
                ch_label = 'unk ch'
                print('could not extract channel label form roi set name')
            
            for radius_ind, radius in enumerate(unique_radius_values):
                current_color = sns.color_palette()[radius_ind]
            #query = {condition: u_value}
            #trials = filterTrials(roi_data.get('epoch_response'), ImagingData, query)
                y = np.mean(on_center_mean_response_interp[:,ch_ind,u_ind,radius_ind,:],axis=0).T
                error = scipy.stats.sem(on_center_mean_response_interp[:,ch_ind,u_ind,radius_ind,:],axis=0).T
                ax1[ch_ind,u_ind].plot(time_vector,y, color = current_color, label = 'r={}'.format(radius))
                ax1[ch_ind, u_ind].legend(loc='upper left')
                ax1[ch_ind,u_ind].fill_between(time_vector, y-error, y+error, color = current_color, alpha = 0.3) #, linestyle='-', color=ImagingData.colors[0])
                ax1[ch_ind,u_ind].set_title('{}, Intensity = {}, {}ms Flash'.format(ch_label,u_value,1000*run_parameters['stim_time']))
                ax1[ch_ind, u_ind].set_ylabel('Response (dF/F)')
                ax1[ch_ind, u_ind].set_xlabel('Time (s)')
                ax1[ch_ind, u_ind].axvspan(run_parameters['pre_time'], run_parameters['pre_time'] + run_parameters['stim_time'], color='gray', alpha=0.2)

def plotAllResponsesByCondition_DS(ImagingDataObjects, ch_names, response_set_name, condition, roi_prefix='rois'):
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
    func_channels_num = len(ch_names)

    for exp_ind, ImagingData in enumerate(ImagingDataObjects):
        fly_metadata = ImagingData.getSubjectMetadata()
        print('fly_metadata: ' + repr(fly_metadata))

        # IMPORTANT: SET TIMING_CHANNEL_IND (visanalysis assumes 0 and will default to first photodiode (ie fly left) of not set!)
        prep = fly_metadata['prep']
    
        if prep == 'fly right optic lobe': #TODO: add condition to check if there are multiple PD channels?  Currently assuming there are two PD recordings
            timing_channel_ind = 1
        elif prep == 'fly left optic lobe':
            timing_channel_ind = 0
        else:
            'could not find photodiode channel based on prep, defaulting to 0'
            timing_channel_ind = 0
        print('photodiode timing_channel_ind: ' + repr(timing_channel_ind))
        ImagingData.timing_channel_ind = timing_channel_ind # IMPORTANT: set timing channel index for photodiode
        for ch_ind, ch_name in enumerate(ch_names):
            # get roi data
            roi_data[exp_ind,ch_ind] = ImagingData.getRoiResponses(ch_name, roi_prefix=roi_prefix)
            # extract mean_response and unique_parameter_values by condition 
            
            unique_parameter_values, mean_response[exp_ind,ch_ind], sem_response[exp_ind,ch_ind], trial_response_by_stimulus[exp_ind,ch_ind] = ImagingData.getTrialAverages(roi_data[exp_ind,ch_ind]['epoch_response'], parameter_key='angle')

            #print('roi_data["epoch_response"].shape: {}'.format(roi_data['epoch_response'].shape))
            #print('roi_data["roi_response"][0].shape: {}'.format(roi_data['roi_response'][0].shape))
            #print('roi_data["time_vector"].shape: {}'.format(roi_data['time_vector'].shape))
            #roi_data['epoch_response'] - 3D array of responses for each trial (roi, trial, time)
            #roi_data['roi_response'] - list of roi responses (roi, time)
            #roi_data['time_vector'] - 1d array of timepoints (one set for all measurements)
            #print('type(roi_data): {}'.format(type(roi_data)))
            #print('roi_data.keys(): {}'.format(roi_data.keys()))
            #print('.time_vector: {}'.format(roi_data.get('time_vector')))
            if ch_ind==0: #only need to do this once for first channel
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

    ds_mean_response_interp = np.empty((n_roi_total, func_channels_num, 2, n_timepoints)) # numpy arrays(roi x channel x direction (PD=0 ND = 1) x time)
    roi_counter = 0
    for exp_ind, ImagingData in enumerate(ImagingDataObjects):
        
        # aggregate PD,ND mean data for each roi

        roi_directions = ImagingData.getRoiParameters(roi_set_name=response_set_name, parameter='direction', roi_prefix='aligned')
        print('roi_directions: ' + repr(roi_directions))


        n_roi = mean_response[exp_ind,0].shape[0] #n_roi in experiment, using first channel as reference
        for roi_ind in range(n_roi): 
            
            p_direction = roi_directions[roi_ind]
            # calculate null direction from preferred
            if p_direction >=  180:
                n_direction = p_direction - 180
            elif p_direction < 180:
                n_direction = p_direction + 180

            p_ind = 0
            n_ind = 1

            for ch_ind, ch_name in enumerate(ch_names):
                for u_ind, up in enumerate(unique_parameter_values):
                    current_direction = up[0]
                    if current_direction == p_direction:
                       
                       # interpolate mean_response and sem_response to common time vector
                        f = interp1d(roi_data[exp_ind,ch_ind]['time_vector'],mean_response[exp_ind,ch_ind][roi_ind,u_ind,:],kind='linear',axis=0,bounds_error = False)
                        ds_mean_response_interp[roi_ind + roi_counter, ch_ind, p_ind,:] = f(time_vector)

                    elif current_direction == n_direction:
                        # interpolate mean_response and sem_response to common time vector
                        f = interp1d(roi_data[exp_ind,ch_ind]['time_vector'],mean_response[exp_ind,ch_ind][roi_ind,u_ind,:],kind='linear',axis=0,bounds_error = False)
                        ds_mean_response_interp[roi_ind + roi_counter, ch_ind, n_ind,:] = f(time_vector)
                        
                    
        roi_counter = roi_counter + n_roi

    print('ds_mean_response_interp.shape (roi, channel, direction (PD, ND), time): ' + repr(ds_mean_response_interp.shape))

    fig_format = '.pdf'
    fig_name_string = 'mean_responses_PD_ND'

    fh1, ax1 = plt.subplots(len(ch_names), 2, figsize=(10, 10*9/16),constrained_layout = True)
    for u_ind, u_value in enumerate(unique_parameter_values):
        for ch_ind, ch_name in enumerate(ch_names):
            if 'ch1' in ch_name:
                ch_label = 'ch1'
            elif 'ch2' in ch_name:
                ch_label = 'ch2'
            else:
                ch_label = 'unk ch'
                print('could not extract channel label form roi set name')
            if u_ind == 0:
                current_direction = 'Preferred'
            elif u_ind == 1:
                current_direction = 'Null'
            #query = {condition: u_value}
            #trials = filterTrials(roi_data.get('epoch_response'), ImagingData, query)
            y = np.mean(ds_mean_response_interp[:,ch_ind,u_ind,:],axis=0).T
            error = scipy.stats.sem(ds_mean_response_interp[:,ch_ind,u_ind,:],axis=0).T
            ax1[ch_ind,u_ind].plot(time_vector,y, color = 'k')
            ax1[ch_ind,u_ind].fill_between(time_vector, y-error, y+error, color = 'c', alpha = 0.2) #, linestyle='-', color=ImagingData.colors[0])
            ax1[ch_ind,u_ind].set_title('{}, {}'.format(ch_label,current_direction))
            ax1[ch_ind, u_ind].set_ylabel('Response (dF/F)')
            ax1[ch_ind, u_ind].set_xlabel('Time (s)')
            ax1[ch_ind, u_ind].axvspan(run_parameters['pre_time'], run_parameters['pre_time'] + run_parameters['stim_time'], color='gray', alpha=0.2)


def plotAllResponsesByConditionComparison(ImagingDataObjects, ch_names, condition, roi_prefix='rois'):
    # plot all roi responses by 
    # unpack roi_data and unique_parameter_values from all ImagingDataObjects

    #ImagingDataObjects is a list of lists where each higher level list is an experiment type
    # that contains ImagingDataObject instances for that experiment type
    #ch_names is a list of roi names (ex ['mask_ch1', 'mask_ch2'])

    roi_data={}
    unique_parameter_values=[]
    mean_response={}
    sem_response={}
    trial_response_by_stimulus={}
    sample_period = []
    n_timepoints = []
    n_roi_total = 0
    for cond_ind, Condition in enumerate(ImagingDataObjects):
        for exp_ind, ImagingData in enumerate(ImagingDataObjects[cond_ind]):
            
            fly_metadata = ImagingData.getSubjectMetadata()
            print('fly_metadata: ' + repr(fly_metadata))

            # IMPORTANT: SET TIMING_CHANNEL_IND (visanalysis assumes 0 and will default to first photodiode (ie fly left) of not set!)
            prep = fly_metadata['prep']
        
            if prep == 'fly right optic lobe': #TODO: add condition to check if there are multiple PD channels?  Currently assuming there are two PD recordings
                timing_channel_ind = 1
            elif prep == 'fly left optic lobe':
                timing_channel_ind = 0
            else:
                'could not find photodiode channel based on prep, defaulting to 0'
                timing_channel_ind = 0
            print('photodiode timing_channel_ind: ' + repr(timing_channel_ind))
            ImagingData.timing_channel_ind = timing_channel_ind # IMPORTANT: set timing channel index for photodiode

            for ch_ind, ch_name in enumerate(ch_names):
                # get roi data
                roi_data[cond_ind,exp_ind,ch_ind] = ImagingData.getRoiResponses(ch_name, roi_prefix=roi_prefix)
                # extract mean_response and unique_parameter_values by condition 
                unique_parameter_values, mean_response[cond_ind,exp_ind,ch_ind], sem_response[cond_ind,exp_ind,ch_ind], trial_response_by_stimulus[cond_ind,exp_ind,ch_ind] = ImagingData.getTrialAverages(roi_data[cond_ind,exp_ind,ch_ind]['epoch_response'], parameter_key='intensity')

                #print('roi_data["epoch_response"].shape: {}'.format(roi_data['epoch_response'].shape))
                #print('roi_data["roi_response"][0].shape: {}'.format(roi_data['roi_response'][0].shape))
                #print('roi_data["time_vector"].shape: {}'.format(roi_data['time_vector'].shape))
                #roi_data['epoch_response'] - 3D array of responses for each trial (roi, trial, time)
                #roi_data['roi_response'] - list of roi responses (roi, time)
                #roi_data['time_vector'] - 1d array of timepoints (one set for all measurements)
                #print('type(roi_data): {}'.format(type(roi_data)))
                #print('roi_data.keys(): {}'.format(roi_data.keys()))
                #print('.time_vector: {}'.format(roi_data.get('time_vector')))
                n_roi = mean_response[cond_ind,exp_ind,ch_ind].shape[0]
                n_roi_total = n_roi_total + n_roi
                n_timepoints.append(len(roi_data[cond_ind,exp_ind,ch_ind]['time_vector']))
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

    mean_response_interp = {}

    for cond_ind, Condition in enumerate(ImagingDataObjects):
        for exp_ind, ImagingData in enumerate(ImagingDataObjects[cond_ind]):
            for ch_ind, ch_name in enumerate(ch_names):
                print('cond_ind: {}'.format(cond_ind))
                print('exp_ind: {}'.format(exp_ind))
                print('ch_ind: {}'.format(ch_ind))
                # interpolate
                #f = interp1d(roi_data[exp_ind,roi_ind]['time_vector'],roi_data[exp_ind,roi_ind]['epoch_response'],kind='linear',axis=2)
                f = interp1d(roi_data[cond_ind,exp_ind,ch_ind]['time_vector'],mean_response[cond_ind,exp_ind,ch_ind][:,:,:],kind='linear',axis=2,bounds_error = False)
                if ch_ind==0:
                    response_interp_temp = np.expand_dims(f(time_vector),axis=-1) # add channel dim
                else:
                    response_interp_temp = np.append(response_interp_temp, np.expand_dims(f(time_vector),axis=-1),axis=-1) # add channel dim, expand along ch axis
                print('response_interp_temp.shape: ' + repr(response_interp_temp.shape))
                #mean_response_interp[:,:,:,ch_ind] = response_interp_temp
            if exp_ind==0:
                    # first assignment
                    mean_response_interp[cond_ind] = response_interp_temp
            else:
                    mean_response_interp[cond_ind] = np.append(mean_response_interp[cond_ind],response_interp_temp, axis=0)

    print('cnt mean_response_interp.shape: ' + repr(mean_response_interp[0].shape))
    print('test mean_response_interp.shape: ' + repr(mean_response_interp[1].shape))
    fh1, ax1 = plt.subplots(len(ch_names), len(unique_parameter_values), figsize=(10, 10*9/16),constrained_layout = True)
    for cond_ind, Condition in enumerate(ImagingDataObjects):
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
                y = np.mean(mean_response_interp[cond_ind][:,u_ind,:,ch_ind],axis=0).T
                print('y.shape: ' + repr(y.shape))
                error = scipy.stats.sem(mean_response_interp[cond_ind][:,u_ind,:,ch_ind],axis=0).T
                print('error.shape: ' + repr(error.shape))
                print('time_vector.shape: ' + repr(time_vector.shape))
                fill_color = ['c','darkorange']
                plt.rc('font', size=16) 
                ax1[ch_ind,u_ind].plot(time_vector,y, color = 'k')
                ax1[ch_ind,u_ind].fill_between(time_vector, y-error, y+error, color=fill_color[cond_ind], alpha = 0.4) #, linestyle='-', color=ImagingData.colors[0])
                ax1[ch_ind,u_ind].set_title('{}, Intensity = {}, {}ms Flash'.format(ch_label,u_value,1000*run_parameters['stim_time']))
                ax1[ch_ind, u_ind].set_ylabel('Response (dF/F)')
                ax1[ch_ind, u_ind].set_xlabel('Time (s)')
                ax1[ch_ind, u_ind].axvspan(run_parameters['pre_time'], run_parameters['pre_time'] + run_parameters['stim_time'], color='gray', alpha=0.2)

def plotF0ByConditionComparison(voxel_mean, color , quiet=True):
    """
    Plot average F0 by condition for multiple groups of ImagingDataObjects.

    Parameters
    ----------


    Returns
    -------
    None.

    """
    
    keys = list(voxel_mean.keys())
    channels = list(dict.fromkeys([j for i in voxel_mean.values() for j in i.keys()]))

    if not quiet:
        print('keys: ' + str(keys))
        print('channels: ' + str(channels))
    for ch in channels:
        y={}
        fig = plt.figure() 
        sns.set_theme(font_scale=2.6, style="whitegrid")
        #plt.title('voxel intensity violin plot ' + str(ch))
        plt.ylabel('voxel intensity')
        # make violin plot of voxel_mean for each channel, key
        for key in keys:
            y[key]=voxel_mean[key][ch]
            print('median ({}), {}: '.format(key,ch) + repr(np.median(y[key])))
        sns.violinplot(y, palette=color)
        plt.show()
        plt.close()


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

def filterPklFiles(data_directory,
                    file_search_string='*.pkl',
                    quiet=False,
                    recursive=False):
    """
    Searches through a directory of visprotocol datafiles and finds datafiles/series that match the search values
    Can search based on any number of fly metadata params or run parameters

    Params
        -data_directory: directory of visprotocol data files to search through

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
    return fileNames

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
