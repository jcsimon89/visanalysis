# -*- coding: utf-8 -*-
"""
Bruker / Prairie View plugin for visanalysis.

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
import os
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import h5py
import skimage.io as io
import functools
import nibabel as nib

from visanalysis import plugin


class BrukerPlugin(plugin.base.BasePlugin):
    def __init__(self):
        super().__init__()
        self.current_series = None
        self.mean_brain = None
        self.current_series_number = 0

    def updateImageSeries(self, data_directory, image_file_name, series_number, channel):
        if series_number != self.current_series_number:  # only re-load if selected new series
            self.current_series_number = series_number
            self.loadImageSeries(data_directory=data_directory,
                                 image_file_name=image_file_name,
                                 channel=channel)

    def getRoiImage(self, data_directory, image_file_name, series_number, channel, z_slice):
        if series_number != self.current_series_number:
            self.current_series_number = series_number
            self.loadImageSeries(data_directory=data_directory,
                                 image_file_name=image_file_name,
                                 channel=channel)
        else:
            pass  # don't need to re-load the entire series

        if self.current_series is None:  # No image file found
            roi_image = []
        else:
            roi_image = self.mean_brain[:, :, int(z_slice)]

        return roi_image

    def getRoiDataFromPath(self, roi_path):
        """
        Compute roi response from roi path objects.

        param:
            roi_path: list of path objects

        *Must first define self.current_series
        """
        mask = self.getRoiMaskFromPath(roi_path)

        roi_response = np.mean(self.current_series[mask, :], axis=0, keepdims=True) - np.min(self.current_series)

        return roi_response

    def getRoiMaskFromPath(self, roi_path):
        """
        Compute roi mask from roi path objects.

        param:
            roi_path: list of path objects

        *Must first define self.current_series
        """
        x_dim, y_dim, z_dim, t_dim = self.current_series.shape

        pixX = np.arange(y_dim)
        pixY = np.arange(x_dim)
        xv, yv = np.meshgrid(pixX, pixY)
        roi_pix = np.vstack((xv.flatten(), yv.flatten())).T

        mask = np.zeros(shape=(x_dim, y_dim, z_dim))

        for path in roi_path:
            z_level = path.z_level
            xy_indices = np.reshape(path.contains_points(roi_pix, radius=0.5), (x_dim, y_dim))
            mask[xy_indices, z_level] = 1

        mask = mask == 1  # convert to boolean for masking

        return mask

    def attachData(self, experiment_file_name, file_path, data_directory):
        for series_number in self.getSeriesNumbers(file_path):
            # # # # Retrieve metadata from files in data directory # # #
            file_basename = 'TSeries-' + experiment_file_name.replace('-', '') + '-' + ('00' + str(series_number))[-3:]
            metadata_filepath = os.path.join(data_directory, file_basename)
            if os.path.exists(metadata_filepath + '.xml'):
                # Photodiode trace
                v_rec_suffix = '_Cycle00001_VoltageRecording_001'
                photodiode_basepath = os.path.join(data_directory, file_basename + v_rec_suffix)
                frame_monitor, time_vector, sample_rate = getPhotodiodeSignal(photodiode_basepath)

                # Metadata & timing information
                response_timing = getAcquisitionTiming(metadata_filepath)
                metadata = getMetaData(metadata_filepath)

                # # # # Attach metadata to epoch run group in data file # # #\
                def find_series(name, obj, sn):
                    target_group_name = 'series_{}'.format(str(sn).zfill(3))
                    if target_group_name in name:
                        return obj

                with h5py.File(file_path, 'r+') as experiment_file:
                    find_partial = functools.partial(find_series, sn=series_number)
                    epoch_run_group = experiment_file.visititems(find_partial)

                    # make sure subgroups exist for stimulus and response timing
                    stimulus_timing_group = epoch_run_group.require_group('stimulus_timing')
                    plugin.base.overwriteDataSet(stimulus_timing_group, 'frame_monitor', frame_monitor)
                    plugin.base.overwriteDataSet(stimulus_timing_group, 'time_vector', time_vector)
                    stimulus_timing_group.attrs['sample_rate'] = sample_rate

                    acquisition_group = epoch_run_group.require_group('acquisition')
                    plugin.base.overwriteDataSet(acquisition_group, 'time_points', response_timing['stack_times'])
                    if 'frame_times' in response_timing:
                        plugin.base.overwriteDataSet(acquisition_group, 'frame_times', response_timing['frame_times'])
                    acquisition_group.attrs['sample_period'] = response_timing['sample_period']
                    for key in metadata:
                        acquisition_group.attrs[key] = metadata[key]

                print('Attached data to series {}'.format(series_number))
            else:
                print('WARNING! Required metadata files not found at {}'.format(metadata_filepath))

    def loadImageSeries(self, data_directory, image_file_name, channel=0):
        metadata_image_dims = self.ImagingDataObject.getAcquisitionMetadata().get('image_dims')  # xyztc
        image_file_path = os.path.join(data_directory, image_file_name)
        if '.tif' in image_file_name:  # tif assumed to be tyx series
            # native axes order is tyx: convert to xyzt, with z dummy axis
            image_series = io.imread(image_file_path)
            image_series = np.swapaxes(image_series, 0, 2)[:, :, np.newaxis, :]  # -> xyzt
            print('Loaded xyt image series {}'.format(image_file_path))
        elif '.nii' in image_file_name:
            nib_brain = np.asanyarray(nib.load(image_file_path).dataobj)
            brain_dims = nib_brain.shape
            if len(brain_dims) == 3:  # xyt
                image_series = nib_brain[:, :, np.newaxis, :]  # -> xyzt
                print('Loaded xyt image series {}'.format(image_file_path))

            elif len(brain_dims) == 4:  # xyzt or xytc
                if brain_dims[-1] == metadata_image_dims[-1]:  # xytc
                    image_series = np.squeeze(nib_brain[:, :, :, channel])[:, :, np.newaxis, :]  # xytc -> xyzt
                    print('Loaded xytc image series {}'.format(image_file_path))
                else:  # xyzt
                    image_series = nib_brain  # xyzt
                    print('Loaded xyzt image series {}'.format(image_file_path))

            elif len(brain_dims) == 5:  # xyztc
                image_series = np.squeeze(nib_brain[:, :, :, :, channel])  # xyzt
                print('Loaded xyzt image series from xyztc {}: channel {}'.format(image_file_path, channel))

            else:
                print('Unrecognized image dimensions')
                image_series = None
        else:
            print('Unrecognized image format. Expects .tif or .nii')

        self.current_series = image_series
        self.mean_brain = np.mean(image_series, axis=3)  # xyz

    # %%
    ###########################################################################
    # Functions for timing and metadata
    #   Accessible outside of the plugin object
    ###########################################################################


def getPhotodiodeSignal(filepath):
    """

    params:
        :filepath: path to photodiode file(s), with no suffix
    """

    metadata = ET.parse(filepath + '.xml')
    root = metadata.getroot()
    rate_node = root.find('Experiment').find('Rate')
    sample_rate = int(rate_node.text)

    active_channels = []
    signal_list = list(root.find('Experiment').find('SignalList'))
    for signal_node in signal_list:
        is_channel_active = signal_node.find('Enabled').text
        channel_name = signal_node.find('Name').text
        if is_channel_active == 'true':
            active_channels.append(channel_name)

    # Load frame tracker signal and pull frame/epoch timing info
    data_frame = pd.read_csv(filepath + '.csv')

    time_vector = data_frame.get('Time(ms)').values / 1e3  # ->sec

    frame_monitor = []  # get responses in all active channels
    for ac in active_channels:
        frame_monitor.append(data_frame.get(' ' + ac).values)
    frame_monitor = np.vstack(frame_monitor)

    return frame_monitor, time_vector, sample_rate


def getMetaData(filepath):
    """

    params:
        filepath: path to photodiode file(s), with no suffix

    returns
        metadata: dict

    """
    metaData = ET.parse(filepath + '.xml')
    root = metaData.getroot()

    metadata = {}
    for child in list(root.find('PVStateShard')):
        if child.get('value') is None:
            for subchild in list(child):
                if subchild.get('value') is None:
                    for subsubchild in list(subchild):
                        new_key = child.get('key') + '_' + subchild.get('index') + subsubchild.get('subindex')
                        new_value = subsubchild.get('value')
                        metadata[new_key] = new_value

                else:
                    new_key = child.get('key') + '_' + subchild.get('index')
                    new_value = subchild.get('value')
                    metadata[new_key] = new_value

        else:
            new_key = child.get('key')
            new_value = child.get('value')
            metadata[new_key] = new_value

    # Get axis dims
    sequences = root.findall('Sequence')
    c_dim = len(sequences[0].findall('Frame')[0].findall('File'))  # number of channels
    x_dim = metadata['pixelsPerLine']
    y_dim = metadata['linesPerFrame']

    if root.find('Sequence').get('type') == 'TSeries Timed Element':  # Plane time series
        t_dim = len(sequences[0].findall('Frame'))
        z_dim = 1
    elif root.find('Sequence').get('type') == 'TSeries ZSeries Element':  # Volume time series
        t_dim = len(sequences)
        z_dim = len(sequences[0].findall('Frame'))
    elif root.find('Sequence').get('type') == 'ZSeries':  # Single Z stack (anatomical)
        t_dim = 1
        z_dim = len(sequences[0].findall('Frame'))
    else:
        print('!Unrecognized series type in PV metadata!')

    metadata['image_dims'] = [int(x_dim), int(y_dim), z_dim, t_dim, c_dim]

    metadata['version'] = root.get('version')
    metadata['date'] = root.get('date')
    metadata['notes'] = root.get('notes')

    return metadata


def getAcquisitionTiming(filepath):
    """
    Imaging acquisition metadata based on the bruker metadata file (xml)

    params:
        filepath: path to photodiode file(s), with no suffix

    returns
        response_timing: dict

    """
    metaData = ET.parse(filepath + '.xml')
    root = metaData.getroot()

    if root.find('Sequence').get('type') == 'TSeries ZSeries Element':
        # volumetric xyz time series
        num_t = len(root.findall('Sequence'))
        num_z = len(root.find('Sequence').findall('Frame'))
        frame_times = np.ndarray(shape=(num_t, num_z), dtype=float)
        frame_times[:] = np.nan
        for t_ind, t_step in enumerate(root.findall('Sequence')):
            for z_ind, z_step in enumerate(t_step.findall('Frame')):
                frame_times[t_ind, z_ind] = z_step.get('relativeTime')

        stack_times = frame_times[:, 0]
        sample_period = np.mean(np.diff(stack_times))

        response_timing = {'frame_times': frame_times,
                           'stack_times': stack_times,
                           'sample_period': sample_period}

    elif root.find('Sequence').get('type') == 'TSeries Timed Element':
        # Single-plane, xy time series
        stack_times = []
        for frame in root.find('Sequence').findall('Frame'):
            frTime = frame.get('relativeTime')
            stack_times.append(float(frTime))

        stack_times = np.array(stack_times)

        sample_period = np.mean(np.diff(stack_times))  # sec
        response_timing = {'stack_times': stack_times,
                           'sample_period': sample_period}

    return response_timing


def get_mark_points_metadata(filepath):
    """
    Parse Bruker / PrairieView markpoints metadata from .xml file.

    params:
        filepath: path to photodiode file(s), with no suffix

    returns
        metadata: dict
    """
    metadata = {}

    root = ET.parse(filepath + '.xml').getroot()
    for key in root.keys():
        metadata[key] = root.get(key)

    point_element = root.find('PVMarkPointElement')
    for key in point_element.keys():
        metadata[key] = point_element.get(key)

    galvo_element = point_element[0]
    for key in galvo_element.keys():
        metadata[key] = galvo_element.get(key)

    points = list(galvo_element)
    for point_ind, point in enumerate(points):
        for key in point.keys():
            metadata['Point_{}_{}'.format(point_ind+1, key)] = point.get(key)

    return metadata
