"""test meanbrain sharpness with different numbers of frames included
"""
import numpy as np
import nibabel as nib
import os
import matplotlib.pyplot as plt


# load example series
path = os.path.join('C:/Users/jcsimon/Documents/Stanford/Data/Bruker/eyesss/JS140_x_JS251/fly_006/func0/imaging/channel_2.nii')
slice = np.asarray(nib.load(path).dataobj[:,:,0,:],dtype='float32') # load first slice
print(slice.shape)
# calculate mean for different length time chunks
chunk_length = [4, 9, 19, 49, 99, 199, 399, 799, 1599, 3199, 4463]
for ind, length in enumerate(chunk_length):
    mean_img = np.mean(slice[:,:,:length], axis=2) # time
    #plot
    fig = plt.subplot(1,1,1)
    fig.imshow(mean_img)
    fig.set_title('frames: {}'.format(length + 1))
    plt.show()
