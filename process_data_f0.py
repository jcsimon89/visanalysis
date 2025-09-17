"""
process data:
1. attach metadata
2. draw roi masks with reduced gui
3. extract roi responses and save to hdf5

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""

import sys
import os
import argparse
import json
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from visanalysis.plugin import base as base_plugin
from visanalysis.analysis import imaging_data
from visanalysis.util import plot_tools
import nibabel as nib
import seaborn as sns

# call structure: python process_data.py --experiment_file_directory "path" --rig "rigID" --series_number "series_number" --run_gui "True/False" --attach_metadata "True/False"

# experiment_file_directory: string that contains the path to the folder that has the .hdf5 file
# rig: rigID (Bruker or AODscope)

#NOTE: all file names and paths are relative to the experiment_file_directory which is the main fly folder ie /fly_001


if __name__ == '__main__':
    ## parse shell arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment_file_directory", nargs="?", help="Folder pointing to hdf5")
    parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
    parser.add_argument("--run_gui", nargs="?", help="True/False")
    parser.add_argument("--attach_metadata", nargs="?", help="True/False")
    parser.add_argument("--show_figs", nargs="?", help="True/False")
    parser.add_argument("--save_figs", nargs="?", help="True/False")
    parser.add_argument("--series_number", nargs="?", help="series numbers for functional data")
    parser.add_argument("--channel_number", nargs="*", help="channel number for each series")
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig
    
    if args.show_figs == 'True':
        show_figs = True
    elif args.show_figs == 'False':
        show_figs = False
    else:
        show_figs = False
        print('not able to interperet show_figs flag, must be "True" or "False", default = False')
    print('show_figs: ' + str(show_figs))

    if args.save_figs == 'True':
        save_figs = True
    elif args.save_figs == 'False':
        save_figs = False
    else:
        save_figs = False
        print('not able to interperet save_figs flag, must be "True" or "False", default = False')
    print('save_figs: ' + str(save_figs))

    if args.run_gui == 'True':
        run_gui = True
    elif args.run_gui == 'False':
        run_gui = False
    else:
        print('not able to interperet run_gui flag, must be "True" or "False"')
    
    if args.attach_metadata == 'True':
        attach_metadata = True
    elif args.attach_metadata == 'False':
        attach_metadata = False
    else:
        print('not able to interperet attach_metadata flag, must be "True" or "False"')

    # series and channel numbers for functional data

    series_num = args.series_number
    channel_num = args.channel_number

    print('series_num = ' + str(series_num))
    print('channel_num = ' + str(channel_num))

    # hardcoded names
    experiment_file_name = 'fly.hdf5' #name of hdf5 file
    json_file_name = 'fly.json' #name of json file
    roi_set_name = 'roi_set_name' #name you're going to save rois as in gui

    # derive image_file_path
    # assumes folder structure from snake_eyesss

    # file_path (complete path to .hdf5 file)
    experiment_file_path = os.path.join(experiment_file_directory, experiment_file_name)

    #TODO extract subject number from hdf5 or folder name? not sure if needed

    print('experiment_file_path: ' + repr(os.path.join(experiment_file_directory, experiment_file_name)))

    ## determine correct visanalysis plugin

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

    # attach metadata to hdf5
    print('attach_metadata: ' + str(attach_metadata))
    if attach_metadata:
        plug.attachData(experiment_file_name, experiment_file_directory)

        print('Attached metadata to {}'.format(experiment_file_name))

    ### loop through experiments (one per optic lobe)

    roi_mask_bool = {} # dict with key = exp and value = bool array of roi masks, shape: roi_index, x, y ,(z)
    image_data = {} # dict with key = exp,ch and value = 4D array of image data, shape: x,y,z,t
    mean_image = {} # dict with key = exp,ch and value = 3D array of mean image data, shape: x,y,z
    voxel_mean = {} # dict with key = ch and value = flat list of voxel values in roi mask (not separated by roi)
    

    for exp_ind, current_series_num in enumerate(series_num): #loop through all series, current_series_num is a list of series to analyze together
        exp='exp' + str(exp_ind)
        current_channel_num = channel_num[exp_ind]

        # define filenames for drawing rois and extracting masks

        roi_channel = int(current_channel_num[0]) # methods expect channel number to be datatype int
        roi_series = int(current_series_num[0]) # use first element of tuple to draw rois for current series, methods expect series number to be datatype int
        
        image_file_name = 'channel_' + str(roi_channel) + '.nii' 
        image_relative_directory = 'func' + str(roi_series) + '/imaging' #folder where .nii is, assumes func_ folder counting starts from 0 which series counter starts from 1
        image_file_directory = os.path.join(experiment_file_directory, image_relative_directory)
        image_file_path = os.path.join(image_file_directory, image_file_name)
        print('drawing rois on series ' + str(roi_series) + 'using channel ' + str(roi_channel))
        print('image_file_path: ' + str(image_file_path))
        
        ##draw roi masks using reduced GUI
        print('run_gui: ' + str(run_gui))
        if run_gui:

            gui_path = str(os.path.join(os.getcwd(),"gui/DataGUI_prog.py"))

            os.system('python ' + gui_path
                    + ' --experiment_file_directory ' + experiment_file_directory
                    + ' --experiment_file_name ' + experiment_file_name
                    + ' --experiment_file_path ' + experiment_file_path
                    + ' --rig ' + rig
                    + ' --series_number ' + str(roi_series) # use first element of tuple to draw rois for current series
                    + ' --image_file_path ' + image_file_path)


        ID = imaging_data.ImagingDataObject(experiment_file_path,
                                            roi_series, # use first element of tuple to extract rois for current series
                                            quiet=False)

        # Retrieve saved mask region responses from data file

        roi_mask_bool[exp] = ID.getRoiMasks(roi_set_name)['roi_mask'] #shape: roi_index, x, y ,(z)
        roi_image = ID.getRoiMasks(roi_set_name)['roi_image'] #shape: x,y,(z)

        print('dimensions of roi_image: ' + str(np.shape(roi_image)))
        print('dimensions of roi_mask_bool[exp]: ' + str(np.shape(roi_mask_bool[exp])))

        # figure out if data is volume or slice
        num_spacial_dim = len(roi_image.shape)

        # add func_spatial_dim and get fly metadata for later
        fly_metadata = ID.getSubjectMetadata()
        print('fly_metadata: ' + repr(fly_metadata))
    
        # convert roi mask to unique numbers for each separate roi (background = 0, first roi =1, second roi =2 etc.)
        # this is how masks need to be formatted for plug.saveRegionResponsesFromMask()
        roi_mask = np.zeros(roi_mask_bool[exp].shape[1:])
        print('dimensions of func data: ' + str(np.shape(roi_mask))) # shape:x,y,(z)

        if num_spacial_dim == 3: #data is a volume
            for roi_ind in range(roi_mask_bool[exp].shape[0]):
                for i in range(roi_mask_bool[exp].shape[1]):
                    for j in range(roi_mask_bool[exp].shape[2]):
                        for k in range(roi_mask_bool[exp].shape[3]):
                            if roi_mask_bool[exp][roi_ind,i,j,k] == True:
                                roi_mask[i,j,k]=roi_ind+1 #since roi_ind starts at 0 and we want to label the first roi with value=1
        else:
            print('data does not appear to be a volume, num_spatial_dim = ' + str(num_spacial_dim))

        n_roi = len(np.unique(roi_mask))-1 # numpy integer
        print('number of rois in roi_mask: ' + str(n_roi))

        ## loop through series and extract mean series image

        for series_ind, current_series in enumerate(current_series_num):
            current_series = int(current_series) # methods expect series number to be datatype int
            sn = 'sn' + str(current_series)
            current_channel = int(current_channel_num[series_ind]) # methods expect channel number to be datatype int
            ch = 'ch' + str(current_channel)

            image_file_name = 'channel_' + str(current_channel) + '.nii' 
            image_relative_directory = 'func' + str(current_series) + '/imaging' #folder where .nii is, assumes func_ folder counting starts from 0 which series counter starts from 1
            image_file_directory = os.path.join(experiment_file_directory, image_relative_directory)
            image_file_path = os.path.join(image_file_directory, image_file_name)
            print('processing series ' + str(current_series))
            print('image_file_path: ' + str(image_file_path))

            # load image data, take mean across time dimension

            image_proxy = nib.load(image_file_path)
            image_data[exp,ch] = np.asarray(image_proxy.dataobj, dtype=np.float32) # shape: x,y,z,t
            mean_image[exp,ch] = np.mean(image_data[exp,ch], axis=3) # shape: x,y,z

            # save mean image file if it doesn't already exist
            mean_image_file_name = 'channel_' + str(current_channel) + '_mean.nii'
            mean_image_file_path = os.path.join(image_file_directory, mean_image_file_name)
            if not os.path.exists(mean_image_file_path):
                mean_image_nifti = nib.Nifti1Image(mean_image[exp,ch], image_proxy.affine, image_proxy.header)
                nib.save(mean_image_nifti, mean_image_file_path)
                print('saved mean image to ' + str(mean_image_file_path))
            else:
                print('mean image file already exists at ' + str(mean_image_file_path))

            # extract flat voxel mean intensity for all rois from mean image
            
            if mean_image.shape == roi_mask.shape: #check that mean image and roi mask have same dimensions
                if ch not in voxel_mean:
                    voxel_mean[ch] = mean_image[exp,ch][roi_mask!=0]
                else:
                    np.append(voxel_mean[ch], mean_image[exp,ch][roi_mask!=0]) 
            else:
                print('mean image and roi mask do not have same dimensions, series ' + str(current_series))
            


            ## 2. plots

            # 1: all roi overlay on mean image

            fig_name_string = 'mask_overlay_all'
            fig_format = '.pdf'
            fig_dir = os.path.join(experiment_file_directory, 'figs')

            for slice in range(roi_mask_bool[exp].shape[3]): #loop through z slices
                im = plot_tools.overlayImage(mean_image[exp,ch][:,:,slice], roi_mask_bool[exp][:,:,:,slice], alpha=0.5, colors=None)
                fig = plt.subplot(1,1,1)
                fig.imshow(im)
                if save_figs:
                    fig_name = fig_name_string + '_rois_{}_slice_{}_'.format(exp,slice)
                    plt.savefig(os.path.join(fig_dir,fig_name + fig_format), dpi=400, transparent=True)
                if show_figs:
                    plt.show()
                plt.close()
    
        # 2: violin plot of voxel intensity
        fig_name_string = 'voxel_intensity_violin'
        fig_format = '.pdf'
        fig_dir = os.path.join(experiment_file_directory, 'figs')
        for ch in voxel_mean:
            fig = plt.figure()
            sns.violinplot(y=voxel_mean[ch])
            plt.title('voxel intensity violin plot ' + str(ch))
            plt.ylabel('voxel intensity')
            if save_figs:
                fig_name = fig_name_string + '_' + str(ch)
                plt.savefig(os.path.join(fig_dir,fig_name + fig_format), dpi=400, transparent=True)
            if show_figs:
                plt.show()
            plt.close()

