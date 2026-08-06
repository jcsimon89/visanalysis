#!/usr/bin/env python3
"""
Data file GUI and ROI drawer.

https://github.com/ClandininLab/visanalysis
mhturner@stanford.edu
"""
import os
import psutil
import sys
import argparse

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import path
import matplotlib.colors as mcolors
from matplotlib.widgets import LassoSelector, EllipseSelector
import matplotlib.cm as cm
from PyQt6.QtWidgets import (QPushButton, QWidget, QLabel, QGridLayout,
                             QApplication, QComboBox, QLineEdit,
                             QSlider, QMessageBox)
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
import numpy as np

from visanalysis.util import plot_tools, h5io


class DataGUI(QWidget):

    def __init__(self):
        super().__init__()
        # parse shell arguments
        parser = argparse.ArgumentParser()
        parser.add_argument("--experiment_file_directory", nargs="?", help="Folder pointing to hdf5")
        parser.add_argument("--experiment_file_name", nargs="?", help="experiment_file_name.hdf5")
        parser.add_argument("--experiment_file_path", nargs="?", help="complete path to .hdf5 file")
        parser.add_argument("--rig", nargs="?", help="Bruker or AODscope")
        parser.add_argument("--series_number", nargs="?", help="integer starting with 1")
        parser.add_argument("--image_file_path", nargs="?", help="full path to .nii image file to be used for roi selection")
        args = parser.parse_args()
        experiment_file_directory = args.experiment_file_directory
        experiment_file_name = args.experiment_file_name
        experiment_file_path = args.experiment_file_path
        rig = args.rig
        series_number = args.series_number
        image_file_path = args.image_file_path
        print('experiment_file_directory = ' + str(experiment_file_directory))
        print('experiment_file_name = ' + str(experiment_file_name))
        print('rig = ' + str(rig))
        self.experiment_file_name = experiment_file_name #name of hdf5 file, previously without .hdf5 but i'm changing that
        self.experiment_file_directory = experiment_file_directory #path to folder that contains .hdf5 file
        self.experiment_file_path = experiment_file_path #complete path to .hdf5 file
        self.rig = rig 
        self.image_file_path = image_file_path
        self.image_file_name = os.path.split(image_file_path)[-1]
        self.image_file_directory = os.path.split(image_file_path)[0]
        self.max_rois = 999
        self.roi_type = 'freehand'
        self.roi_radius = None
        self.existing_roi_set_paths = {}
        self.current_roi_index = 0
        self.current_z_slice = 0
        self.current_channel = 1  # index
        self.series_number = series_number
        self.roi_response = []
        self.roi_mask = []
        self.roi_path = []
        self.roi_image = None
        self.fano_image = None
        self.show_fano_overlay = False
        self.roi_path_list = []

        self.blank_image = np.zeros((1, 1))

        self.colors = [mcolors.to_rgb(x) for x in list(mcolors.XKCD_COLORS)[:self.max_rois]]

        self.fano_cmap = plt.get_cmap('viridis').copy()
        self.fano_cmap.set_bad(alpha=0)  # drawn-over ROI pixels (NaN) show as transparent

        self.img_data_min, self.img_data_max = 0.0, 1.0
        self.img_vmin, self.img_vmax = 0.0, 1.0
        self.fano_data_min, self.fano_data_max = 0.0, 1.0
        self.fano_vmin, self.fano_vmax = 0.0, 1.0

        self.initializeDataAnalysis()

        self.initUI()

        self.plugin.updateImagingDataObject(experiment_file_directory, experiment_file_name, series_number)

        self.updateExistingRoiSetList()
        self.selectImageDataFile()

    def initUI(self):
        self.grid = QGridLayout(self)

        self.file_info_grid = QGridLayout()
        self.file_info_grid.setSpacing(3)
        self.grid.addLayout(self.file_info_grid, 0, 0, 1, 2)

        self.roi_control_grid = QGridLayout()
        self.roi_control_grid.setSpacing(3)
        self.grid.addLayout(self.roi_control_grid, 0, 2)

        self.plot_grid = QGridLayout()
        self.plot_grid.setSpacing(3)
        self.grid.addLayout(self.plot_grid, 1, 2)

        # Label with current expt file
        self.currentExperimentLabel = QLabel('')
        self.currentExperimentLabel.setText(self.experiment_file_name)
        self.file_info_grid.addWidget(self.currentExperimentLabel, 0, 0)

        self.experiment_file_path_display = QLabel('')
        self.experiment_file_path_display.setText('..' + self.experiment_file_path[-24:])
        self.experiment_file_path_display.setFont(QtGui.QFont('SansSerif', 8))
        self.file_info_grid.addWidget(self.experiment_file_path_display, 1, 0)

        # File name display
        self.currentImageFileNameLabel = QLabel('')
        self.file_info_grid.addWidget(self.currentImageFileNameLabel, 2, 0)

        # # # # Roi control # # # # # # # # (0, 2)
        # ROI type drop-down
        self.RoiTypeComboBox = QComboBox(self)
        self.RoiTypeComboBox.addItem("freehand")
        radii = [1, 2, 3, 4, 6, 8]
        for radius in radii:
            self.RoiTypeComboBox.addItem("circle:"+str(radius))
        self.RoiTypeComboBox.activated.connect(self.selectRoiType)
        self.roi_control_grid.addWidget(self.RoiTypeComboBox, 0, 0)

        # Clear all ROIs button
        self.clearROIsButton = QPushButton("Clear ROIs", self)
        self.clearROIsButton.clicked.connect(self.clearRois)
        self.roi_control_grid.addWidget(self.clearROIsButton, 0, 2)

        # Response display type dropdown
        self.RoiResponseTypeComboBox = QComboBox(self)

        self.RoiResponseTypeComboBox.addItem("RawTrace")
        self.RoiResponseTypeComboBox.addItem("TrialAverage")
        self.RoiResponseTypeComboBox.addItem("TrialResponses")
        self.RoiResponseTypeComboBox.addItem("TrialAverageDFF")
        self.RoiResponseTypeComboBox.addItem("TrialAverageBright")
        self.roi_control_grid.addWidget(self.RoiResponseTypeComboBox, 2, 2)

        # ROIset file name line edit box
        self.defaultRoiSetName = "roi_set_name"
        self.le_roiSetName = QLineEdit(self.defaultRoiSetName)
        self.roi_control_grid.addWidget(self.le_roiSetName, 1, 1)

        # Save ROIs button
        self.saveROIsButton = QPushButton("Save ROIs", self)
        self.saveROIsButton.clicked.connect(self.saveRois)
        self.roi_control_grid.addWidget(self.saveROIsButton, 1, 0)

        # Load ROI set combobox
        self.loadROIsComboBox = QComboBox(self)
        self.loadROIsComboBox.addItem("(load existing ROI set)")
        self.loadROIsComboBox.activated.connect(self.selectedExistingRoiSet)
        self.roi_control_grid.addWidget(self.loadROIsComboBox, 1, 2)
        self.updateExistingRoiSetList()

        # Delete current roi button
        self.deleteROIButton = QPushButton("Delete ROI", self)
        self.deleteROIButton.clicked.connect(self.deleteRoi)
        self.roi_control_grid.addWidget(self.deleteROIButton, 2, 0)

        # Current roi slider
        self.roiSlider = QSlider(QtCore.Qt.Orientation.Horizontal, self)
        self.roiSlider.setMinimum(0)
        self.roiSlider.setMaximum(self.max_rois)
        self.roiSlider.valueChanged.connect(self.sliderUpdated)
        self.roi_control_grid.addWidget(self.roiSlider, 2, 1, 1, 1)

        ctx = plt.rc_context({'xtick.major.size': 1,
                              'axes.spines.top': False,
                              'axes.spines.right': False,
                              'xtick.labelsize': 'xx-small',
                              'ytick.labelsize': 'xx-small',
                              'xtick.major.size': 1.0,
                              'ytick.major.size': 1.0,
                              'xtick.major.pad': 1.0,
                              'ytick.major.pad': 1.0})
        with ctx:
            self.responseFig = plt.figure(frameon=False, layout='constrained')
            self.responsePlot = self.responseFig.add_subplot(111)
            self.responseCanvas = FigureCanvas(self.responseFig)
        self.responseCanvas.draw_idle()
        self.plot_grid.addWidget(self.responseCanvas, 0, 0, 1, 2)

        # # # # Image canvas # # # # # # # # (1, 2)
        self.roi_fig = plt.figure()
        self.roi_ax = self.roi_fig.add_subplot(111)
        self.roi_canvas = FigureCanvas(self.roi_fig)
        self.toolbar = NavigationToolbar(self.roi_canvas, self)
        self.roi_ax.set_aspect('equal')
        self.roi_ax.set_axis_off()
        self.plot_grid.addWidget(self.toolbar, 1, 0)
        self.plot_grid.addWidget(self.roi_canvas, 2, 0)

        # Fano overlay colorbar: create once so its reserved space never resizes roi_ax
        # on toggle -- shown/hidden afterward, never removed/re-added.
        overlay_im_placeholder = self.roi_ax.imshow(self.blank_image, cmap=self.fano_cmap, alpha=0.5)
        self.overlay_cbar = self.roi_fig.colorbar(overlay_im_placeholder, ax=self.roi_ax, fraction=0.046, pad=0.04)
        self.overlay_cbar.ax.set_visible(False)

        self.roi_canvas.mpl_connect('scroll_event', self.onScrollZoom)

        # # # # Fano factor canvas (variance / mean over time) # # # # (1, 2)
        self.fano_fig = plt.figure()
        self.fano_ax = self.fano_fig.add_subplot(111)
        self.fano_canvas = FigureCanvas(self.fano_fig)
        self.fano_ax.set_aspect('equal')
        self.fano_ax.set_axis_off()
        self.fano_ax.set_title('Fano factor', fontsize=8)
        self.fano_im = self.fano_ax.imshow(self.blank_image, cmap=self.fano_cmap)
        self.fano_cbar = self.fano_fig.colorbar(self.fano_im, ax=self.fano_ax, fraction=0.046, pad=0.04)
        self.plot_grid.addWidget(self.fano_canvas, 2, 1)

        self.plot_grid.setRowStretch(0, 1)
        self.plot_grid.setRowStretch(1, 3)
        self.plot_grid.setRowStretch(2, 3)
        self.plot_grid.setColumnStretch(0, 1)
        self.plot_grid.setColumnStretch(1, 1)

        # Toggle fano factor overlay on roi image (below the roi image, not the fano panel)
        self.fanoOverlayButton = QPushButton("Show Fano Overlay", self)
        self.fanoOverlayButton.setCheckable(True)
        self.fanoOverlayButton.toggled.connect(self.toggleFanoOverlay)
        self.plot_grid.addWidget(self.fanoOverlayButton, 3, 0)

        # Contrast (min/max) sliders, one pair per canvas
        self.imgMinSlider, self.imgMaxSlider, img_contrast_grid = self._buildContrastSliders('Image', self.imgContrastChanged)
        self.plot_grid.addLayout(img_contrast_grid, 4, 0)

        self.fanoMinSlider, self.fanoMaxSlider, fano_contrast_grid = self._buildContrastSliders('Fano', self.fanoContrastChanged)
        self.plot_grid.addLayout(fano_contrast_grid, 4, 1)

        # Current z slice slider
        self.zSlider = QSlider(QtCore.Qt.Orientation.Horizontal, self)
        self.zSlider.setMinimum(0)
        self.zSlider.setMaximum(50)
        self.zSlider.setValue(0)
        self.zSlider.valueChanged.connect(self.zSliderUpdated)
        self.plot_grid.addWidget(self.zSlider, 5, 0, 1, 2)

        self.roi_fig.tight_layout()

        self.setWindowTitle('Visanalysis')
        self.setGeometry(200, 200, 1200, 600)
        self.show()

    def _buildContrastSliders(self, label, on_change):
        """Build a labeled min/max contrast slider pair (0-1000, mapped to a data range elsewhere).

        min slider stacked above max slider; each slider's left end is its low
        value and right end is its high value.
        """
        grid = QGridLayout()
        min_slider = QSlider(QtCore.Qt.Orientation.Horizontal, self)
        min_slider.setRange(0, 1000)
        min_slider.setValue(0)
        min_slider.valueChanged.connect(on_change)
        max_slider = QSlider(QtCore.Qt.Orientation.Horizontal, self)
        max_slider.setRange(0, 1000)
        max_slider.setValue(1000)
        max_slider.valueChanged.connect(on_change)
        grid.addWidget(QLabel('{} min'.format(label)), 0, 0)
        grid.addWidget(min_slider, 0, 1)
        grid.addWidget(QLabel('{} max'.format(label)), 1, 0)
        grid.addWidget(max_slider, 1, 1)
        return min_slider, max_slider, grid

    def updateExistingRoiSetList(self):
        if self.experiment_file_name is not None:
            file_path = self.experiment_file_path
            self.existing_roi_set_paths = self.plugin.getRoiSetPaths(file_path)  # dictionary of name: full path
            self.loadROIsComboBox.clear()
            for r_path in self.existing_roi_set_paths:
                self.loadROIsComboBox.addItem(r_path)

            self.show()

    def selectedExistingRoiSet(self):
        file_path = self.experiment_file_path
        roi_set_key = self.loadROIsComboBox.currentText()
        roi_set_path = self.existing_roi_set_paths[roi_set_key]

        _, _, self.roi_path, self.roi_mask = self.plugin.loadRoiSet(file_path, roi_set_path)

        if self.series_number is not None:
            self.roi_response = []
            for new_path in self.roi_path:
                new_roi_resp = self.plugin.getRoiDataFromPath(roi_path=new_path)
                self.roi_response.append(new_roi_resp)

            # update slider to show most recently drawn roi response
            self.current_roi_index = len(self.roi_response)-1
            self.roiSlider.setValue(self.current_roi_index)

            # Update figures
            self.redrawRoiTraces()

    def initializeDataAnalysis(self):
        file_path = self.experiment_file_path
        data_type = self.rig
        # Load plugin based on Rig name in hdf5 file
        if data_type == 'Bruker':
            from visanalysis.plugin import bruker
            self.plugin = bruker.BrukerPlugin()
        elif data_type == 'AODscope':
            from visanalysis.plugin import aodscope
            self.plugin = aodscope.AodScopePlugin()
        else:
            from visanalysis.plugin import base
            self.plugin = base.BasePlugin()

        self.plugin.parent_gui = self

        # # # TEST # # #
        memory_usage = psutil.Process(os.getpid()).memory_info().rss*10**-9
        print('Current memory usage: {:.2f}GB'.format(memory_usage))
        sys.stdout.flush()
        # # # TEST # # #

    def selectImageDataFile(self):
        file_path = self.experiment_file_path
        print('User selected image file at {}'.format(self.image_file_path))
        h5io.attachImageFileName(file_path, self.series_number, self.image_file_name)
        print('Attached image_file_name {} to series {}'.format(self.image_file_name, self.series_number))
        print('Data directory is {}'.format(self.experiment_file_path))

        self.currentImageFileNameLabel.setText(self.image_file_name)

        # show roi image
        if self.series_number is not None:
            if self.experiment_file_path is not None:  # user has selected a raw data directory
                self.plugin.updateImageSeries(data_directory=self.image_file_directory,
                                              image_file_name=self.image_file_name,
                                              series_number=self.series_number,
                                              channel=self.current_channel)
                self.roi_image = self.plugin.mean_brain
                self.fano_image = self.computeFanoImage()
                self.resetContrastRanges()
                self.zSlider.setValue(0)
                self.zSlider.setMaximum(self.roi_image.shape[2]-1)
                self.redrawRoiTraces(reset_zoom=True)
            else:
                print('Select a data directory before drawing rois')

# %% # # # # # # # # ROI SELECTOR WIDGET # # # # # # # # # # # # # # # # # # #

    def refreshLassoWidget(self, keep_paths=False, reset_zoom=False):
        preserve_zoom = self.roi_image is not None and not reset_zoom
        if preserve_zoom:
            prev_xlim = self.roi_ax.get_xlim()
            prev_ylim = self.roi_ax.get_ylim()

        self.roi_ax.clear()
        init_lasso = False
        if self.roi_image is not None:
            image_slice = self.roi_image[:, :, self.current_z_slice]
            if len(self.roi_mask) > 0:
                newImage = plot_tools.overlayImage(image_slice, self.roi_mask, 0.5, self.colors,
                                                   z=self.current_z_slice, vmin=self.img_vmin, vmax=self.img_vmax)
                self.roi_ax.imshow(newImage)
            else:
                self.roi_ax.imshow(image_slice, cmap=cm.gray, vmin=self.img_vmin, vmax=self.img_vmax)
            init_lasso = True

            if self.show_fano_overlay and self.fano_image is not None:
                fano_slice = self.getFanoSliceForDisplay(self.current_z_slice)
                # Re-point the persistent colorbar at a fresh overlay image rather than
                # adding/removing a colorbar, which would resize roi_ax on every toggle.
                overlay_im = self.roi_ax.imshow(fano_slice, cmap=self.fano_cmap, alpha=0.5,
                                                vmin=self.fano_vmin, vmax=self.fano_vmax)
                self.overlay_cbar.update_normal(overlay_im)
                self.overlay_cbar.ax.set_visible(True)
            else:
                self.overlay_cbar.ax.set_visible(False)

            if preserve_zoom:  # imshow() resets view limits to fit the new image -- restore prior zoom/pan
                self.roi_ax.set_xlim(prev_xlim)
                self.roi_ax.set_ylim(prev_ylim)
        else:
            self.roi_ax.imshow(self.blank_image)
            self.overlay_cbar.ax.set_visible(False)
        self.roi_ax.set_axis_off()

        self.roi_canvas.draw()

        if not keep_paths:
            self.roi_path_list = []

        if init_lasso:
            if self.roi_type == 'circle':
                self.lasso_1 = EllipseSelector(self.roi_ax, onselect=self.newEllipse, button=1)
            elif self.roi_type == 'freehand':
                self.lasso_1 = LassoSelector(self.roi_ax, onselect=self.newFreehand, button=1)
                self.lasso_2 = LassoSelector(self.roi_ax, onselect=self.appendFreehand, button=3)
            else:
                print('Warning ROI type not recognized. Choose circle or freehand')

        self.refreshFanoWidget()

    def onScrollZoom(self, event):
        """Zoom the roi image in/out around the cursor. Doesn't touch data coordinates,
        so it can't affect where a lasso/ellipse ROI actually lands."""
        if event.inaxes != self.roi_ax or self.roi_image is None:
            return
        if event.xdata is None or event.ydata is None:
            return
        zoom_factor = 0.8 if event.button == 'up' else 1.25
        xlim = self.roi_ax.get_xlim()
        ylim = self.roi_ax.get_ylim()
        self.roi_ax.set_xlim([event.xdata + (x - event.xdata) * zoom_factor for x in xlim])
        self.roi_ax.set_ylim([event.ydata + (y - event.ydata) * zoom_factor for y in ylim])
        self.roi_canvas.draw_idle()

    def refreshFanoWidget(self):
        if self.fano_image is not None:
            fano_slice = self.getFanoSliceForDisplay(self.current_z_slice)
            self.fano_im.set_data(fano_slice)
            self.fano_im.set_extent((-0.5, fano_slice.shape[1] - 0.5, fano_slice.shape[0] - 0.5, -0.5))
            self.fano_ax.set_xlim(-0.5, fano_slice.shape[1] - 0.5)
            self.fano_ax.set_ylim(fano_slice.shape[0] - 0.5, -0.5)
            self.fano_im.set_clim(vmin=self.fano_vmin, vmax=self.fano_vmax)
        else:
            self.fano_im.set_data(self.blank_image)
        self.fano_canvas.draw()

    def getFanoSliceForDisplay(self, z_slice):
        """Fano factor image for one z slice, with already-drawn ROI pixels excluded (NaN)."""
        fano_slice = self.fano_image[:, :, z_slice].copy()
        if len(self.roi_mask) > 0:
            drawn = np.any([mask[:, :, z_slice] for mask in self.roi_mask], axis=0)
            fano_slice[drawn] = np.nan
        return fano_slice

    def toggleFanoOverlay(self, checked):
        self.show_fano_overlay = checked
        self.refreshLassoWidget(keep_paths=True)

    def resetContrastRanges(self):
        """Auto-scale contrast sliders to the newly loaded image/fano data."""
        self.img_data_min = float(np.nanmin(self.roi_image))
        self.img_data_max = float(np.nanmax(self.roi_image))
        self.img_vmin, self.img_vmax = self.img_data_min, self.img_data_max
        self.imgMinSlider.blockSignals(True)
        self.imgMaxSlider.blockSignals(True)
        self.imgMinSlider.setValue(0)
        self.imgMaxSlider.setValue(1000)
        self.imgMinSlider.blockSignals(False)
        self.imgMaxSlider.blockSignals(False)

        finite_fano = self.fano_image[np.isfinite(self.fano_image)] if self.fano_image is not None else np.array([])
        if finite_fano.size > 0:
            self.fano_data_min = float(np.min(finite_fano))
            self.fano_data_max = float(np.max(finite_fano))
        else:
            self.fano_data_min, self.fano_data_max = 0.0, 1.0
        self.fano_vmin, self.fano_vmax = self.fano_data_min, self.fano_data_max
        self.fanoMinSlider.blockSignals(True)
        self.fanoMaxSlider.blockSignals(True)
        self.fanoMinSlider.setValue(0)
        self.fanoMaxSlider.setValue(1000)
        self.fanoMinSlider.blockSignals(False)
        self.fanoMaxSlider.blockSignals(False)

    def _sliderToValue(self, slider, data_min, data_max):
        frac = slider.value() / slider.maximum()
        return data_min + frac * (data_max - data_min)

    def imgContrastChanged(self):
        if self.imgMinSlider.value() >= self.imgMaxSlider.value():
            if self.sender() is self.imgMinSlider:
                self.imgMaxSlider.setValue(self.imgMinSlider.value() + 1)
            else:
                self.imgMinSlider.setValue(self.imgMaxSlider.value() - 1)
        self.img_vmin = self._sliderToValue(self.imgMinSlider, self.img_data_min, self.img_data_max)
        self.img_vmax = self._sliderToValue(self.imgMaxSlider, self.img_data_min, self.img_data_max)
        self.refreshLassoWidget(keep_paths=True)

    def fanoContrastChanged(self):
        if self.fanoMinSlider.value() >= self.fanoMaxSlider.value():
            if self.sender() is self.fanoMinSlider:
                self.fanoMaxSlider.setValue(self.fanoMinSlider.value() + 1)
            else:
                self.fanoMinSlider.setValue(self.fanoMaxSlider.value() - 1)
        self.fano_vmin = self._sliderToValue(self.fanoMinSlider, self.fano_data_min, self.fano_data_max)
        self.fano_vmax = self._sliderToValue(self.fanoMaxSlider, self.fano_data_min, self.fano_data_max)
        self.refreshLassoWidget(keep_paths=True)

    def newFreehand(self, verts):
        new_roi_path = path.Path(verts)
        new_roi_path.z_level = self.zSlider.value()
        new_roi_path.channel = self.current_channel
        self.updateRoiSelection([new_roi_path])

    def appendFreehand(self, verts):
        print('Appending rois, hit Enter/Return to finish')
        new_roi_path = path.Path(verts)
        new_roi_path.z_level = self.zSlider.value()
        new_roi_path.channel = self.current_channel
        self.roi_path_list.append(new_roi_path)

    def keyPressEvent(self, event):
        if type(event) == QtGui.QKeyEvent:
            if np.any([event.key() == QtCore.Qt.Key.Key_Return, event.key() == QtCore.Qt.Key.Key_Enter]):
                if len(self.roi_path_list) > 0:
                    event.accept()
                    self.updateRoiSelection(self.roi_path_list)
                else:
                    event.ignore()
            else:
                event.ignore()
        else:
            event.ignore()

    def newEllipse(self, pos1, pos2, definedRadius=None):
        x1 = np.round(pos1.xdata)
        x2 = np.round(pos2.xdata)
        y1 = np.round(pos1.ydata)
        y2 = np.round(pos2.ydata)

        radiusX = np.sqrt((x1 - x2)**2)/2
        radiusY = np.sqrt((y1 - y2)**2)/2
        if self.roi_radius is not None:
            radiusX = self.roi_radius

        center = (np.round((x1 + x2)/2), np.round((y1 + y2)/2))
        new_roi_path = path.Path.circle(center=center, radius=radiusX)
        new_roi_path.z_level = self.zSlider.value()
        new_roi_path.channel = self.current_channel
        self.updateRoiSelection([new_roi_path])

    def updateRoiSelection(self, new_roi_path):
        mask = self.plugin.getRoiMaskFromPath(new_roi_path)
        new_roi_resp = self.plugin.getRoiDataFromPath(roi_path=new_roi_path)
        if mask.sum() == 0:
            print('No pixels in the roi you just drew')
            return
        # update list of roi data
        self.roi_mask.append(mask)
        self.roi_path.append(new_roi_path)  # list of lists of paths
        self.roi_response.append(new_roi_resp)
        # update slider to show most recently drawn roi response
        self.current_roi_index = len(self.roi_response)-1
        self.roiSlider.setValue(self.current_roi_index)

        # Update figures
        self.redrawRoiTraces()

    def sliderUpdated(self):
        self.current_roi_index = self.roiSlider.value()
        self.redrawRoiTraces()

    def zSliderUpdated(self):
        self.current_z_slice = self.zSlider.value()
        if self.roi_image is not None:
            self.refreshLassoWidget(keep_paths=True)

    def redrawRoiTraces(self, reset_zoom=False):
        self.clearRoiArtists()
        if self.current_roi_index < len(self.roi_response):
            current_raw_trace = np.squeeze(self.roi_response[self.current_roi_index])
            fxn_name = self.RoiResponseTypeComboBox.currentText()
            display_trace = getattr(self.plugin, 'getRoiResponse_{}'.format(fxn_name))([current_raw_trace])
            self.responsePlot.plot(display_trace, color=self.colors[self.current_roi_index], linewidth=1, alpha=0.5)
            self.responsePlot.set_xlim([0, len(display_trace)])
            y_min = np.nanmin(display_trace)
            y_max = np.nanmax(display_trace)
            self.responsePlot.set_ylim([y_min, y_max])
        self.responseCanvas.draw()

        self.refreshLassoWidget(keep_paths=False, reset_zoom=reset_zoom)

# %% # # # # # # # # LOADING / SAVING / COMPUTING ROIS # # # # # # # # # # # # # # # # # # #

    def loadRois(self, roi_set_path):
        file_path = self.experiment_file_path
        self.roi_response, self.roi_image, self.roi_path, self.roi_mask = self.plugin.loadRoiSet(file_path, roi_set_path)
        self.zSlider.setValue(0)
        self.zSlider.setMaximum(self.roi_image.shape[2]-1)


    def saveRois(self):
        file_path = self.experiment_file_path
        roi_set_name = self.le_roiSetName.text()
        if roi_set_name in h5io.getAvailableRoiSetNames(file_path, self.series_number):
            buttonReply = QMessageBox.question(self,
                                               'Overwrite roi set',
                                               "Are you sure you want to overwrite roi set: {}?".format(roi_set_name),
                                               QMessageBox.StandardButton.Yes |
                                               QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
            if buttonReply == QMessageBox.StandardButton.Yes:
                self.plugin.saveRoiSetMask(file_path,
                                       series_number=self.series_number,
                                       roi_set_name=roi_set_name,
                                       roi_mask=self.roi_mask,
                                       #roi_response=self.roi_response,
                                       roi_image=self.roi_image,
                                       roi_path=self.roi_path)
                print('Saved roi set {} to series {}'.format(roi_set_name, self.series_number))
                self.updateExistingRoiSetList()
            else:
                print('Overwrite aborted - pick a unique roi set name')
        else:
            self.plugin.saveRoiSetMask(file_path,
                                   series_number=self.series_number,
                                   roi_set_name=roi_set_name,
                                   roi_mask=self.roi_mask,
                                   #roi_response=self.roi_response,
                                   roi_image=self.roi_image,
                                   roi_path=self.roi_path)
            print('Saved roi set {} to series {}'.format(roi_set_name, self.series_number))
            self.updateExistingRoiSetList()

    def deleteRoi(self):
        if self.current_roi_index < len(self.roi_response):
            self.roi_mask.pop(self.current_roi_index)
            self.roi_response.pop(self.current_roi_index)
            self.roi_path.pop(self.current_roi_index)
            self.roiSlider.setValue(self.current_roi_index-1)
            self.redrawRoiTraces()

    def clearRois(self):
        self.roi_mask = []
        self.roi_response = []
        self.roi_path = []
        self.roi_image = None
        self.fano_image = None
        self.clearRoiArtists()
        self.redrawRoiTraces()
        self.roi_ax.clear()

    def clearRoiArtists(self):
        for artist in self.responsePlot.lines + self.responsePlot.collections:
            artist.remove()

    def selectRoiType(self):
        self.roi_type = self.RoiTypeComboBox.currentText().split(':')[0]
        if 'circle' in self.RoiTypeComboBox.currentText():
            self.roi_radius = int(self.RoiTypeComboBox.currentText().split(':')[1])
        else:
            self.roi_radius = None
        self.redrawRoiTraces()

    def computeFanoImage(self):
        """Fano factor (variance / mean) at each pixel, computed over the time axis of the raw series."""
        current_series = getattr(self.plugin, 'current_series', None)
        if current_series is None:
            return None
        mean_image = np.mean(current_series, axis=3)
        var_image = np.var(current_series, axis=3)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(mean_image > 0, var_image / mean_image, np.nan)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = DataGUI()
    sys.exit(app.exec())
