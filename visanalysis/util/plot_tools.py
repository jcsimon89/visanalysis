"""
Assorted plotting utilities.

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
import numpy as np
import matplotlib.pyplot as plt


def addLine(ax, x, y, line_name='', color='k', linestyle='-', marker='None'):
    ax.plot(x, y, linestyle=linestyle, marker=marker,
            linewidth=1, label=line_name, color=color)


def addErrorBars(ax, xdata, ydata, line_name='',
                 stat='sem',
                 mode='sticks',
                 color='k'):

    if len(xdata.shape) == 2 and len(ydata.shape) == 2:  # x and y error
        err_x = _calcError(xdata, stat)
        err_y = _calcError(ydata, stat)
        mean_x = np.mean(xdata, axis=0)
        mean_y = np.mean(ydata, axis=0)

        _addYError(ax, mean_x, mean_y, err_y, line_name, mode, stat, color)
        _addXError(ax, mean_x, mean_y, err_x, line_name, mode, stat, color)

    elif len(xdata.shape) == 1 and len(ydata.shape) == 2:  # y error
        err_y = _calcError(ydata, stat)
        mean_y = np.mean(ydata, axis=0)

        _addYError(ax, xdata, mean_y, err_y, line_name, mode, stat, color)
    elif len(xdata.shape) == 2 and len(ydata.shape) == 1:  # x error
        err_x = _calcError(xdata, stat)
        mean_x = np.mean(xdata, axis=0)

        _addXError(ax, mean_x, ydata, err_x, line_name, mode, stat, color)
    else:
        raise Exception('no population data to compute errors')


def addScaleBars(axis, dT, dF, T_value=-0.1, F_value=-0.4):
    axis.plot(T_value * np.ones((2)), np.array([F_value, F_value + dF]), 'k-', alpha=0.9)
    axis.plot(np.array([T_value, dT + T_value]), F_value * np.ones((2)), 'k-', alpha=0.9)


def _addXError(ax, x, y, err_x, line_name, mode, stat, color):
    xx = ax.plot([x - err_x, x + err_x], [y, y], linestyle='--', marker=None, linewidth=1, label=line_name + '_errX')


def _addYError(ax, x, y, err_y, line_name, mode, stat, color):
    yp = ax.plot(x, y - err_y, linestyle='--', marker=None,
                 linewidth=1, label=line_name + '_errY_plus', color=color)
    ym = ax.plot(x, y + err_y, linestyle='--', marker=None,
                 linewidth=1, label=line_name + '_errY_minus', color=color)
    ym[0].tag = mode
    yp[0].tag = 'hide'


def _calcError(data, stat):
    if stat == 'sem':
        err = np.std(data, axis=0) / np.sqrt(data.shape[0])
    elif stat == 'std':
        err = np.std(data, axis=0)
    return err


# tools for images:
def addImageScaleBar(ax, image, scale_bar_length, microns_per_pixel, location):
    dim_x = image.shape[1]
    dim_y = image.shape[0]
    dx = scale_bar_length / microns_per_pixel  # pixels
    if location[0] == 'l':
        start_y = 0.9 * dim_y
    elif location[0] == 'u':
        start_y = 0.1 * dim_y

    if location[1] == 'l':
        start_x = 0.1 * dim_x
        end_x = start_x + dx
    elif location[1] == 'r':
        start_x = 0.9 * dim_x
        end_x = start_x - dx
    ax.plot([start_x, end_x], [start_y, start_y], 'w')


def overlayImage(im, mask, alpha, colors=None, z=0, vmin=None, vmax=None):
    # image = [x,y,rgb]
    # mask can be 4d with slices
    #mask = [rois,x,y,(z)]
    # mask should be binary (true for mask, false for background)
    # vmin/vmax: contrast range to normalize im by (defaults to [0, im.max()], as before)

    if vmin is None:
        vmin = 0
    if vmax is None:
        vmax = np.max(im)
    im = np.clip((im - vmin) / (vmax - vmin), 0, 1) # normalize image to contrast range
    if len(im.shape) < 3:
        imRGB = np.tile(im[..., np.newaxis], 3) # add rgb vals as 3rd dim
    else:
        imRGB = im

    overlayComponent = 0
    origImageComponent = 0
    if len(mask[0].shape) == 2: # mask[0] is first roi
        compositeMask = np.tile(mask[0][..., np.newaxis], 3) # if mask[0] is 2d, add rgb vals in 3rd dim
    else:
        compositeMask = np.tile(mask[0][:, :, z, np.newaxis], 3) # if mask[0] is 3d, take slice z and then add rgb vals in 3rd dim
    for ind, currentRoi in enumerate(mask):
        if len(mask[0].shape) == 2: 
            maskRGB = np.tile(currentRoi[..., np.newaxis], 3) # mask with rgb in 3rd dim for specific roi
        else:
            maskRGB = np.tile(currentRoi[:, :, z, np.newaxis], 3) # mask with rgb in 3rd dim for specific roi
        if colors is None:
            newColor = (1, 0, 0)
        else:
            newColor = colors[ind]
        
        #compositeMask = [x,y,rgb]
        #maskRGB =  [x,y,rgb] for specific roi

        compositeMask = compositeMask + maskRGB #add each maskRGB to composite mask as you loop through rois
        overlayComponent += alpha * np.array(newColor) * maskRGB # make semitransparent mask overlay for specific roi
        origImageComponent += (1 - alpha) * maskRGB * imRGB # image component in the area where the overlay will happen

    untouched = (compositeMask == False) * imRGB # part of image where no overlay is happening

    im_out = untouched + overlayComponent + origImageComponent # prevents rest of image and underlay from being dimmed during multiple overlays
    im_out = (im_out * 255).astype(np.uint8)
    return im_out


def cleanAxes(ax):
    ax.set_axis_off()
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])


def plotBinnedEpochResponse(ImagingDataObject, roi_data, bin_frequency,
                            roi_inds=None,
                            stim_timing=None,
                            dff='pre',
                            color='k', fill_color=None, fill_alpha=0.3,
                            title=None, ylabel=None, xlabel='Time (s)',
                            figsize=None, ax=None):
    """
    Plot mean ± SEM response across epochs using binned time vectors.

    Calls ImagingDataObject.getBinnedEpochAverage() to bin the epoch response
    matrix onto a common time grid at the specified frequency, then plots
    the mean trace with SEM shading for each requested ROI.

    Params:
        ImagingDataObject: ImagingDataObject instance (provides getBinnedEpochAverage)
        roi_data: dict as returned by getRoiResponses (must contain
            'epoch_response' and 'time_vector_by_epoch')
        bin_frequency: float, Hz. Temporal frequency for the output time grid.
        roi_inds: list of int or None. Which ROIs to plot. None = all ROIs.
        stim_timing: dict with 'pre_time' and 'stim_time' (sec) to shade the
            stimulus window, or None to skip.
        color: line color (any matplotlib color spec)
        fill_color: SEM fill color. None defaults to same as line color.
        fill_alpha: float, transparency of SEM fill (0-1)
        title: str or None, figure title
        ylabel: str, y-axis label
        xlabel: str, x-axis label
        figsize: tuple or None, figure size (width, height)
        ax: matplotlib Axes or array of Axes. If None, a new figure is created.

    Returns:
        fh: figure handle (None if ax was provided)
        axes: array of Axes used
        bin_centers: 1d array, bin center times (sec)
        mean_response: ndarray, shape = (n_rois, n_bins)
        sem_response: ndarray, shape = (n_rois, n_bins)
    """
    # Auto-set ylabel based on dff mode if not explicitly provided
    if ylabel is None:
        if dff == 'none':
            ylabel = 'F'
        else:
            ylabel = r'$\Delta F/F_0$'

    # Compute binned average
    bin_centers, mean_response, sem_response, _ = ImagingDataObject.getBinnedEpochAverage(
        roi_data['epoch_response'],
        roi_data['time_vector_by_epoch'],
        bin_frequency
    )

    n_rois = mean_response.shape[0]
    if roi_inds is None:
        roi_inds = list(range(n_rois))

    if fill_color is None:
        fill_color = color

    # Create figure if no axes provided
    fh = None
    if ax is None:
        if figsize is None:
            figsize = (8, 3 * len(roi_inds))
        fh, ax = plt.subplots(len(roi_inds), 1, figsize=figsize,
                              constrained_layout=True, squeeze=False)
        ax = ax[:, 0]  # flatten to 1D array of axes

    if not hasattr(ax, '__len__'):
        ax = [ax]

    for plot_ind, roi_ind in enumerate(roi_inds):
        cur_ax = ax[plot_ind]
        y = mean_response[roi_ind, :]
        err = sem_response[roi_ind, :]

        cur_ax.plot(bin_centers, y, color=color, linewidth=1.5)
        cur_ax.fill_between(bin_centers, y - err, y + err,
                            color=fill_color, alpha=fill_alpha)

        # Stimulus timing overlay
        if stim_timing is not None:
            pre = stim_timing.get('pre_time', 0)
            stim = stim_timing.get('stim_time', 0)
            cur_ax.axvspan(pre, pre + stim, color='gray', alpha=0.15)

        cur_ax.set_ylabel(ylabel)
        cur_ax.set_xlabel(xlabel)
        if len(roi_inds) > 1:
            cur_ax.set_title('ROI {}'.format(roi_ind))

    if title is not None:
        if fh is not None:
            fh.suptitle(title)
        else:
            ax[0].set_title(title)

    return fh, ax, bin_centers, mean_response, sem_response

