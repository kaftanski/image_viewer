import codecs
import os
from glob import glob
from math import floor, copysign
from typing import Union, List, Sequence

import SimpleITK as sitk
import matplotlib as mpl
import numpy as np
from PyQt5.QtWidgets import QAction
from matplotlib.backend_bases import cursors
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.figure import Figure

try:
    from matplotlib.backend_bases import MouseButton
except ImportError as e:
    # for compatibility with older Matplotlib versions (i.e. 2.1.1), where the MouseButton enum is not yet defined
    import enum

    class MouseButton(enum.IntEnum):
        LEFT = 1
        MIDDLE = 2
        RIGHT = 3

from GUI.viewer_utils import add_mask_to_image, get_aspect_ratio_for_plane, compatible_metadata, index_compatibility, \
    get_3d_plane_index, ImageMask, ImageMarker

mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['image.cmap'] = 'gray'


SLICE_ORIENTATION_XY = 2
SLICE_ORIENTATION_XZ = 1
SLICE_ORIENTATION_YZ = 0

SLICE_ORIENTATION = {'xy': 2, 'xz': 1, 'yz': 0}


class IMIImageViewer(QtWidgets.QMainWindow):
    """ A PyQt window with matplotlib-based viewer of a given image """
    def __init__(self, image: sitk.Image, title: str = '', mask: ImageMask = None):
        """ Constructor with intitial image to be shown

        @param image: SimpleITK Image to be shown
        @param title: optionally specify window title
        @param mask: optional initial mask to be shown over the image
        """
        super(IMIImageViewer, self).__init__()

        self.image_viewer = ImageViewerWidget(image, parent=self)
        if mask is not None:
            self.image_viewer.add_mask(mask)

        self.setCentralWidget(self.image_viewer)
        self.init_menu_bar()
        self.setWindowTitle(title if title != '' else 'IMI Image Viewer')
        self.show()

    def get_markers_for_region_growing(self) -> List[List[int]]:
        """ Get the user input markers in a list of lists of pixel coordinates

        @return: the coordinates of the markers in image coordinates
        """
        return list(list(int(coord) for coord in m.pixel_position) for m in self.image_viewer.markers)

    def init_menu_bar(self):
        """ Creates and populates the PyQt menuBar with viewer control actions
        """
        menubar = self.menuBar()

        file_menu = menubar.addMenu('&File')

        exit_action = QAction('Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close)

        open_action = QAction('&Open Other Image', self)
        open_action.setShortcut('Ctrl+O')
        open_action.setStatusTip('Replace the current image with another image from a file.')
        open_action.triggered.connect(self.load_image)

        open_mask_action = QAction('Load reference &mask', self)
        open_mask_action.setShortcut('Ctrl+M')
        open_mask_action.setStatusTip('Load reference mask of the image')
        open_mask_action.triggered.connect(self.load_mask)

        save_view_action = QAction('&Save Current View', self)
        save_view_action.setShortcut('Ctrl+S')
        save_view_action.setStatusTip('Save the currently displayed image slice as a PNG.')
        save_view_action.triggered.connect(self.save_current_view)

        file_menu.addActions([open_action, open_mask_action, save_view_action, exit_action])

        view_menu = menubar.addMenu('&View')
        orientation_submenu = view_menu.addMenu('Image &Orientation')
        xy_orientation_action = QAction('XY plane', self)
        xy_orientation_action.setStatusTip('Set orientation to XY (axial) plane')
        xy_orientation_action.triggered.connect(lambda: self.image_viewer.change_orientation(SLICE_ORIENTATION_XY))

        xz_orientation_action = QAction('XZ plane', self)
        xz_orientation_action.setStatusTip('Set orientation to XZ (coronal) plane')
        xz_orientation_action.triggered.connect(lambda: self.image_viewer.change_orientation(SLICE_ORIENTATION_XZ))

        yz_orientation_action = QAction('YZ plane', self)
        yz_orientation_action.setStatusTip('Set orientation to YZ (sagittal) plane')
        yz_orientation_action.triggered.connect(lambda: self.image_viewer.change_orientation(SLICE_ORIENTATION_YZ))

        orientation_submenu.addActions([xy_orientation_action, xz_orientation_action, yz_orientation_action])
        view_menu.addSeparator()

        clear_markers_action = QAction('&Clear all markers', self)
        clear_markers_action.setStatusTip('Clear all markers from all images')
        clear_markers_action.triggered.connect(lambda: self.image_viewer.remove_markers(only_last=False))

        remove_last_marker_action = QAction('Remove &last marker', self)
        remove_last_marker_action.setStatusTip('Remove the last placed marker')
        remove_last_marker_action.setShortcut('Ctrl+Z')
        remove_last_marker_action.triggered.connect(lambda: self.image_viewer.remove_markers(only_last=True))

        # marker_submenu = view_menu.addMenu('&Marker')
        view_menu.addActions([clear_markers_action, remove_last_marker_action])
        view_menu.addSeparator()

        reset_wl_action = QAction('&Reset Window / Level', self)
        reset_wl_action.setStatusTip('Reset window-levelling to span the entire image range')
        reset_wl_action.triggered.connect(self.image_viewer.reset_window_level)

        view_menu.addAction(reset_wl_action)

        help_menu = menubar.addMenu('&Help')
        controls_action = help_menu.addAction('&Controls...')
        controls_action.triggered.connect(self.show_controls)

    def show_controls(self):
        """ Open an information dialog detailing the viewer's controls
        """
        message_box = QtWidgets.QMessageBox(self)
        message_box.setWindowTitle('Controls')
        message_box.setTextFormat(QtCore.Qt.RichText)
        with codecs.open((glob('controls.html') + glob('*/controls.html'))[0], 'r') as ctrls_html:
            message_box.setText(ctrls_html.read())

        message_box.exec_()

    def file_dialog(self, caption: str, filter_prompt: str, mode: str = 'save') -> str:
        """ Shows an open or save file dialog and returns the picked path as a string

        @param caption: the string to be displayed as window title
        @param filter_prompt: the file filter prompt (i.e. the file extension, type etc.)
        @param mode: either 'save' or 'open'; determines if an open- or save-file-dialog is shown
        @return the user's chosen file path
        """
        if mode == 'save':
            file_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, caption=caption,
                                                                 directory=os.getcwd(),
                                                                 filter=filter_prompt)
        elif mode == 'open':
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, caption=caption, directory=os.getcwd(),
                                                                 filter=filter_prompt)
        else:
            raise ValueError(
                '\'{}\' is not a valid mode for file_dialog. Choose \'save\' or \'open\' instead.'.format(mode))

        return file_path

    def load_image(self):
        """ Choose a file to open and set it as the new image of the viewer.
        Currently only nifti images are suggested but all image types that ITK can read should work.
        """
        file_path = self.file_dialog(caption='Load Image', filter_prompt='Image (*.nii.gz)', mode='open')

        if file_path != '':
            print('Loading image ' + file_path)

            image = sitk.ReadImage(file_path)
            self.image_viewer.set_image(image)

    def load_mask(self):
        """ Choose a file to open and set it as the red mask over the image.
        Currently only nifti images are suggested but all image types that ITK can read should work.

        It is possible to load any image as mask, because all voxels != 0 are interpreted as object and all
        0-value voxels will be background.
        """
        file_path = self.file_dialog(caption='Load Image Mask', filter_prompt='Image (*.nii.gz)', mode='open')

        if file_path != '':
            print('Loading mask ' + file_path)

            image = sitk.ReadImage(file_path)
            mask = ImageMask(image, color='r')
            self.image_viewer.add_mask(mask)

    def save_current_view(self):
        """ Choose a file to save the current view of the viewer canvas as a .png file.
        """
        file_path = self.file_dialog('Save File', 'Image (*.png, *.PNG)', mode='save')

        if file_path != '':
            if file_path[-4:].lower() != '.png':
                file_path += '.png'

            print('Saving view to {}'.format(file_path))
            self.image_viewer.canvas.print_figure(filename=file_path, dpi=200, bbox_inches=0)


class ImageViewerWidget(QtWidgets.QWidget):
    """ PyQt Widget holding a matplotlib canvas to visualize and navigate through 3D image data
    """
    def __init__(self, image_data: sitk.Image = None, parent: QtWidgets.QWidget = None):
        """ Constructs the viewer and initializes all controls.

        @param image_data: the SimpleITK image to visualize
        @param QtWidgets.QWidget parent: the widget's parent if used in a full PyQt GUI
        """
        super(ImageViewerWidget, self).__init__(parent)
        if image_data is None:
            # create an empty 3d image
            image_data = sitk.GetImageFromArray(np.zeros((1, 1, 1)))

        # init attributes
        self.image = image_data  # the SimpleITK image in the viewer
        self.image_array = None  # the image as np.array
        self.current_slice = 0

        self.im_plot = None  # storing the mpl plot

        self.current_window = 0
        self.current_level = 0

        self.greyval_range = [0, 0]

        self.masks = {}  # dict with mask as key and its plot as value
        self.markers = []  # list of markers and their pixel coordinates
        self.orientation = SLICE_ORIENTATION_XY  # what dimensions are seen as the plane (values like in VTK)

        # init the MPL canvas
        initial_canvas_size_x = 5  # size in pixels. corresponds to a size of 5x5 inches with 100 dpi
        initial_canvas_size_y = 5
        self.max_resolution = 500
        self.dpi = self.max_resolution // initial_canvas_size_x

        self.canvas = FigureCanvas(Figure(figsize=(initial_canvas_size_x, initial_canvas_size_y),
                                          dpi=self.dpi, facecolor='black'))

        self.ax = self.canvas.figure.subplots()
        self.ax.axis('off')
        self.canvas.figure.subplots_adjust(0, 0, 1, 1)

        self.marker_plot = self.ax.scatter([], [], c=ImageMarker.STANDARD_COLOR)

        self.canvas.setParent(self)

        # add a hidden mpl toolbar for panning & zooming functionality
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.hide()

        # define a slider for image slicing
        self.slice_slider = QtWidgets.QSlider(orientation=QtCore.Qt.Vertical, parent=self)
        self.slice_slider.valueChanged.connect(self.update_slice)

        # define a text field to display image coordinates and one for the current zoom
        bottom_bar_layout = QtWidgets.QHBoxLayout()
        self.pixel_info_label = PixelInfoQLabel(parent=self)
        bottom_bar_layout.addWidget(self.pixel_info_label)

        # layout the QtWidget
        main_layout = QtWidgets.QVBoxLayout()
        image_with_slider_layout = QtWidgets.QHBoxLayout()
        image_with_slider_layout.addWidget(self.canvas)
        image_with_slider_layout.addWidget(self.slice_slider)
        main_layout.addLayout(image_with_slider_layout)
        main_layout.addLayout(bottom_bar_layout)
        self.setLayout(main_layout)

        # set the focus to be on the canvas, so that keyevents get handled by mpl_connect
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()

        # prevent the slider from gaining focus
        self.slice_slider.setFocusPolicy(QtCore.Qt.NoFocus)

        # initialize the showing image
        self.init_image()

        # define the interactor for user inputs
        self.interactor_style = ImageViewerInteractor(self)

        # show the widget
        self.show()

    def init_image(self):
        """ Initial call after the image has been set or reset.

        Constructs the image array and plots the middle slice of the image with window/level spanning the whole
        grey value range of the image.
        """
        # DON'T swap axes, in order to keep the fastest dimension of the np array the first
        # -> z-slices are the most performant
        self.image_array = sitk.GetArrayFromImage(self.image)

        # take the middle slice as initial slice shown
        self.current_slice = self.get_slice_dim_size() // 2

        # redraw if there is already a plot
        if self.im_plot is not None:
            self.im_plot.remove()

        self.im_plot = self.ax.imshow(self.image_array[index_compatibility(get_3d_plane_index(self.current_slice, self.orientation))],
                                      aspect=get_aspect_ratio_for_plane(self.image.GetSpacing(), self.orientation, self.image.GetSize()))

        # set the window/level to span the whole greyvalue range of the image
        min_grey = self.image_array.min()
        max_grey = self.image_array.max()
        self.greyval_range = [min_grey, max_grey]
        self.reset_window_level()

        # init the coordinate label
        initial_coords = [0, 0]
        initial_coords.insert(self.orientation, self.current_slice)
        self.show_pixel_info(initial_coords)

        # init the slider with the correct range
        self.adapt_slider()

    def change_orientation(self, orientation: Union[int, str]):
        """ Change the slicing dimension of the viewer.
        Orientation values are expected as in SLICE_ORIENTATION

        @param orientation:
            either the number corresponding to the slicing dimension or a string of the new plane:
                                'xy' or 2 -> the XY plane is shown;
                                'xz' or 1 -> the XZ plane is shown;
                                'yz' or 0 -> the YZ plane is shown;
        """
        # don't do anything if orientation stays the same
        if orientation == self.orientation:
            return

        if orientation not in SLICE_ORIENTATION.keys() and orientation not in SLICE_ORIENTATION.values():
            print('Cannot change viewer orientation to {}. Use one of {} or {} respectively.'.
                  format(orientation, SLICE_ORIENTATION.keys(), SLICE_ORIENTATION.values()))

        # convert the plane string into the number of the slice dimension
        if isinstance(orientation, str):
            orientation = SLICE_ORIENTATION[orientation]

        self.orientation = orientation

        # make sure the current slice is valid in the new orientation
        if self.current_slice >= self.get_slice_dim_size():
            self.current_slice = self.get_slice_dim_size() - 1

        self.redraw_slice()
        self.adapt_slider()

    def set_window_level(self, window: Union[int, float], level: Union[int, float]):
        """ Window levelling operation on the image.

        @param window: width of the window
        @param level: window center
        """
        if window < 0:
            # set a minimum level (could be even smaller in case of normalized images with values in [0,1])
            window = 1e-5

        self.current_window = window
        self.current_level = level

        half_window = self.current_window / 2
        self.im_plot.set_clim(vmin=self.current_level - half_window, vmax=self.current_level + half_window)
        self.canvas.draw()  # optimizable: only redraw the image, not markers or masks!
        # print('Window: {}, Level: {}'.format(self.current_window, self.current_level))

    def reset_window_level(self):
        """ Set the window/level to span the entire grey value range of the image.
        """
        # reset to cover the entire image range
        window = self.greyval_range[1] - self.greyval_range[0]
        level = self.greyval_range[0] + window / 2
        self.set_window_level(window, level)

    def update_slice(self, new_slice: int):
        """ Sets the current slice and redraws the image

        @param new_slice: the slice index
        """
        # if the index is not in range, change it to the min or max index
        lower = 0
        upper = self.get_slice_dim_size()
        if not lower <= new_slice < upper:
            new_slice = lower if new_slice < lower else upper - 1

        self.current_slice = new_slice
        self.redraw_slice()

    def redraw_slice(self):
        """ Draws the current slice of the image with mask and markers.
        This resets the current zoom factor but not the window/level.
        """
        # here, blitting could increase performance

        img_slice = self.image_array[
            index_compatibility(get_3d_plane_index(slice_index=self.current_slice, orientation=self.orientation))]

        # set extent to be according to dimensions so that nothing gets squished
        y_max, x_max = img_slice.shape
        self.ax.set_xlim([0, x_max])
        self.ax.set_ylim([0, y_max])

        half_window = self.current_window / 2

        # remove the image plot and redraw. this is a little less performant than
        # setting data on the current image plot, but this way the extent is updated correctly
        self.im_plot.remove()
        self.im_plot = self.ax.imshow(img_slice, vmin=self.current_level - half_window, vmax=self.current_level + half_window,
                                      aspect=get_aspect_ratio_for_plane(self.image.GetSpacing(), self.orientation, self.image.GetSize()))

        # update the coordinate text
        self.pixel_info_label.set_coordinate(self.orientation, self.current_slice)
        self.show_pixel_info(self.pixel_info_label.coords)

        # draw masks and markers
        self.update_masks()
        self.scatter_markers()

        # update the canvas (making sure this is only called once per update for performance)
        self.canvas.draw()

    def get_slice_dim_size(self) -> int:
        """ Get the size of the slicing dimension (found in self.orientation)

        @return: the size of the current slicing dimension
        """
        return self.image.GetSize()[self.orientation]

    def adapt_slider(self):
        """ Call this, if the image slice range or the current slice changed.
        The slider value and range will be adapted to these values
        """
        self.slice_slider.setRange(0, self.get_slice_dim_size() - 1)
        self.slice_slider.setValue(self.current_slice)

    def move_slice(self, delta: int):
        """ Change the current slice by delta steps and then redraw.

        @param delta: the number of steps by which to change the current slice (positive or negative)
        """
        # moves the slider, which in turn changes the slice and redraws the image
        self.slice_slider.setValue(self.current_slice + delta)

    def add_marker(self, position: Sequence[Union[float, int]]):
        """ Adds a marker to the given pixel position, saves and displays it.

        @param position: continuous 3d pixel position of the new marker
        """
        assert len(position) == 3
        marker = ImageMarker(position)
        self.markers.append(marker)

        # only redraw markers, if the new marker is in the current slice
        if round(position[self.orientation]) == self.current_slice:
            self.scatter_markers()
            self.canvas.draw()

    def remove_markers(self, only_last: bool = False):
        """ Remove all (or only the last) markers from the image.
        This does not only apply to visible markers, but markers in all slices.

        @param only_last: if true, only the last marker will be removed; otherwise all markers get removed
        """
        if len(self.markers) == 0:
            return

        if only_last:
            removed_marker = self.markers.pop()
            print('removed marker at {}'.format(removed_marker.pixel_position))
        else:
            self.markers.clear()
            print('removed all markers')

        # redraw markers
        self.scatter_markers()
        self.canvas.draw()

    def scatter_markers(self):
        """ Draw all markers in the current slice on the image.
        Attention: this does not redraw the whole canvas, self.canvas.draw() has to be called separately.
        This way, multiple changes can be drawn at once, for instance after a new slice is shown.
        """
        # determine which markers are found in the current slice
        markers_in_slice = []
        for m in self.markers:
            pixel_pos = np.array(m.pixel_position)
            # choose a marker if the rounded index is the same as the current slice
            if round(pixel_pos[self.orientation]) == self.current_slice:
                # only take the plane coordinates
                markers_in_slice.append(pixel_pos[np.arange(3) != self.orientation])

        if len(markers_in_slice) != 0:
            # combine marker coordinates into one array and and show it in a scatter plot
            np_markers = np.stack(markers_in_slice, axis=1)
            self.marker_plot.set_offsets(np_markers.swapaxes(0, 1))
        else:
            # if no markers are found, clear remove the current plot and replace it with an empty one
            self.marker_plot.remove()
            self.marker_plot = self.ax.scatter([], [], c=ImageMarker.STANDARD_COLOR, s=ImageMarker.STANDARD_SIZE)

    def add_mask(self, mask: ImageMask):
        """ Show a mask overlay on top of the current image. The mask has to have the same size, spacing and origin.

        @param mask: the mask to overlay (has to have the same properties as the image)
        """
        if not compatible_metadata(self.image, mask.mask_image):
            return

        # save the mask in the masks dict
        self.masks[mask] = None

        # draw the image
        self.update_masks()
        self.canvas.draw()

    def update_masks(self):
        """ Draw the masks in the current slice as an overlay on top of the image.

        Attention: this does not redraw the whole canvas, self.canvas.draw() has to be called separately.
        This way, multiple changes can be drawn at once, for instance after a new slice is shown.
        """
        # remove old plots
        self.clear_mask_plots()

        # draw each mask (no changes to canvas before canvas.draw() is called)
        for m in self.masks.keys():
            self.masks[m] = add_mask_to_image(self.ax, m.get_slice(self.current_slice, self.orientation),
                                              aspect=get_aspect_ratio_for_plane(m.get_spacing(), self.orientation, self.image.GetSize()),
                                              alpha=m.alpha, color=m.color)

    def clear_mask_plots(self):
        """ Remove all masks from the canvas.
        """
        for m in self.masks.keys():
            try:
                # remove the existing plot
                self.masks[m].remove()
            except AttributeError:
                # if the mask was not there in the previous slice it is None here
                pass

    def set_image(self, image: sitk.Image):
        """ Set the image and show it. This resets all markers and masks previously added to the viewer.

        @param image: the new image to show
        """
        self.image = image

        # remove old masks
        self.clear_mask_plots()
        self.masks.clear()

        # remove old markers
        self.markers.clear()
        self.scatter_markers()

        # draw image and re-init everything else
        self.init_image()

    def show_pixel_info(self, pixel_coords: Sequence[Union[float, int]]):
        """ Show the given coordinates and the corresponding image intensity in a label below the image.
        For instance, this is called whenever whenever the user moves the mouse over the image.

        @param pixel_coords: the (continuous) coordinatates to be displayed in the label
        """
        if pixel_coords is None:
            # remove the label text
            self.pixel_info_label.clear()
        else:
            # take the discrete index
            x, y, z = [int(ind) for ind in pixel_coords]
            try:
                self.pixel_info_label.set_values(x, y, z, self.image.GetPixel(x, y, z))
            except (IndexError, RuntimeError, TypeError):
                # gets thrown if x, y, z are out of bounds of the image (this can happen on the edges of the figure)
                return


class PixelInfoQLabel(QtWidgets.QLabel):
    """ Label to show image coordinates and intensity """
    def __init__(self, parent: QtWidgets.QWidget = None):
        """ Initialises the label (with coords=[0,0,0])

        @param parent: the QWidget this element belongs to
        """
        super(PixelInfoQLabel, self).__init__(text='', parent=parent)

        self.coords = [0, 0, 0]
        self.intensity = 0

    def set_values(self, x: Union[float, int], y: Union[float, int], z: Union[float, int], intensity: Union[float, int]):
        """ Updates the label text with given coordinates and intensity

        @param x: x coordinate
        @param y: y coordinate
        @param z: z coordinate
        @param intensity: the image intensity
        """
        self.coords = [x, y, z]
        self.intensity = intensity

        self.update_text()

    def set_coordinate(self, dim, coord):
        """ Set a single coordinate in the given dimension

        @param dim: the dimension the coord belongs to
        @param coord: the index to set
        """
        self.coords[dim] = coord
        self.update_text()

    def set_intensity(self, intensity):
        """ Set the intensity text

        @param intensity: the intensity value to set
        """
        self.intensity = intensity
        self.update_text()

    def update_text(self):
        """ Sets the label text according to this object's attributes
        """
        self.setText('({}, {}, {}): {}'.format(*self.coords, self.intensity))


class ImageViewerInteractor:
    """ Handler for all user inputs on the given ImageViewer """
    def __init__(self, image_viewer: ImageViewerWidget, verbose: bool = False):
        """ Init all connections to events, the mpl canvas can produce

        @param image_viewer: the image viewer to observe
        @param verbose: set True, if you want a console output when actions are triggered
        """
        self.iv = image_viewer
        self.verbose = verbose

        # event observing
        self.iv.canvas.mpl_connect('button_press_event', self.handle_mouse_button_down)
        self.iv.canvas.mpl_connect('button_release_event', self.handle_mouse_button_up)
        self.iv.canvas.mpl_connect('motion_notify_event', self.on_mouse_motion)
        self.iv.canvas.mpl_connect('scroll_event', self.on_mousewheel_event)
        self.iv.canvas.mpl_connect('key_press_event', self.on_key_press)
        # self.iv.canvas.mpl_connect('resize_event', self.on_resize)  # unfinished

        self.last_mouse_position_in_figure = [0, 0]
        self.window_level_start_position = None  # is set on right mouse down in order to determine the new w/l
        self.mouse_wheel_step_residual = 0  # accumulating mousewheel movements, if they dont add up to a whole step

    def on_key_press(self, event: mpl.backend_bases.KeyEvent):
        """ Handles all keys pressed by the user.
        Key -> Functionality:
            up -> move slice forward by 1
            down -> move slice backward by 1
            x -> change orientation to yz plane
            y -> change orientation to xz plane
            z -> change orientation to xy plane
            pageup -> make a page up step through slices (+ 10 by default)
            pagedown -> make a page down step through slices (- 10 by default)

        @param event: the mpl.backend_bases.KeyEvent to handle
        """
        if self.verbose:
            print(event.key)

        if event.key == 'up':
            self.iv.move_slice(1)
        elif event.key == 'down':
            self.iv.move_slice(-1)
        elif event.key == 'x':
            self.iv.change_orientation(SLICE_ORIENTATION['yz'])
        elif event.key == 'y':
            self.iv.change_orientation(SLICE_ORIENTATION['xz'])
        elif event.key == 'z':
            self.iv.change_orientation(SLICE_ORIENTATION['xy'])
        elif event.key == 'r':
            self.iv.reset_window_level()
        elif event.key == 'pageup':
            # emulate the page-up behaviour of the QSlider (as this is never focused)
            self.iv.move_slice(self.iv.slice_slider.pageStep())
        elif event.key == 'pagedown':
            # emulate the page-down behaviour of the QSlider (as this is never focused)
            self.iv.move_slice(- self.iv.slice_slider.pageStep())

    def handle_mouse_button_down(self, event: mpl.backend_bases.MouseEvent):
        """ Handles mouse button down events by distributing event to
        self.on_left_button_down or self.on_right_button_down respectively.

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        if self.verbose:
            print(event.button)

        if event.button == MouseButton.LEFT:
            self.on_left_button_down(event)
        elif event.button == MouseButton.RIGHT:
            self.on_right_button_down(event)
        elif event.button == MouseButton.MIDDLE:
            self.on_middle_button_down(event)

    def handle_mouse_button_up(self, event: mpl.backend_bases.MouseEvent):
        """ Handles mouse button up events by distributing event to
        self.on_left_button_up or self.on_right_button_up respectively.

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        if self.verbose:
            print(event.button, 'up')

        if event.button == MouseButton.LEFT:
            self.on_left_button_up(event)
        elif event.button == MouseButton.RIGHT:
            self.on_right_button_up(event)
        elif event.button == MouseButton.MIDDLE:
            self.on_middle_button_up(event)

    def on_left_button_down(self, event: mpl.backend_bases.MouseEvent):
        """ Left mouse button down event handler:
             if click was inside the figure, a marker is placed at the mouse position

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        if event.inaxes and self.iv.toolbar.mode == '':
            # click inside of figure
            event_data_position = [event.xdata, event.ydata]
            event_data_position.insert(self.iv.orientation, self.iv.current_slice)
            print('adding marker at {}'.format(event_data_position))
            self.iv.add_marker(event_data_position)

    def on_left_button_up(self, event: mpl.backend_bases.MouseEvent):
        """ Left mouse button up event handler:
             no functionality.

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        pass

    def on_middle_button_down(self, event: mpl.backend_bases.MouseEvent):
        """ Left mouse button down event handler:
             toggles zoom/pan mode from the mpl NavigationToolbar and deactivates marker/window-levelling functionality

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        # start zoom/pan mode
        self.iv.toolbar.pan()

        # update cursor manually (otherwise this is only done after the mouse is first moved)
        if self.iv.toolbar.mode == "pan/zoom":
            self.iv.toolbar.set_cursor(cursors.MOVE)  # grabbing hand
        else:
            self.iv.toolbar.set_cursor(cursors.POINTER)  # normal cursor

    def on_right_button_down(self, event: mpl.backend_bases.MouseEvent):
        """ Left mouse button down event handler:
             starts the window levelling.
             Moving the mouse while the right mouse button is pressed does window levelling.
             Movement in x-direction modifies the window, y-direction for level.

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        # if panning / zooming is not activated
        if self.iv.toolbar.mode == '':
            # start window levelling
            self.window_level_start_position = [event.x, event.y]

    def on_right_button_up(self, event: mpl.backend_bases.MouseEvent):
        """ Left mouse button down event handler:
             end of window levelling

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        self.window_level_start_position = None

    def on_middle_button_up(self, event: mpl.backend_bases.MouseEvent):
        """ Left mouse button down event handler:
             no functionality

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        pass

    def on_mouse_motion(self, event: mpl.backend_bases.MouseEvent):
        """ Handles mouse movement events:
             if mouse is over the image, the coordinates and intensity under the cursor are displayed.
             if right mouse button is pressed, window/levelling is carried out
             (see ImageViewerInteractor.on_right_button_down for details)

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        if event.inaxes:
            # update mouse cursor position
            self.last_mouse_position_in_figure = [event.xdata, event.ydata]
            position = [event.xdata, event.ydata]
            position.insert(self.iv.orientation, self.iv.current_slice)
            self.iv.show_pixel_info(position)

        if self.window_level_start_position is not None:
            # window levelling

            # move window by x movement
            # move level by y movement
            dx = event.x - self.window_level_start_position[0]
            dy = event.y - self.window_level_start_position[1]

            if event.x is not None:
                self.window_level_start_position[0] = event.x
            if event.y is not None:
                self.window_level_start_position[1] = event.y

            gain = 5  # gain to amplify or decrease the change
            new_window = self.iv.current_window + dx * gain
            new_level = self.iv.current_level + dy * gain/2

            self.iv.set_window_level(new_window, new_level)

    def on_mousewheel_event(self, event: mpl.backend_bases.MouseEvent):
        """ Handles mouse scrolling events:
             scrolling up moves slice forward, scrolling down moves slice backward.
             scrolling with ctrl key pressed zooms in or out with the cursor position as zoom center.

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        steps = round(event.step)

        # in case of high resolution touchpads/mice: accumulate values under the step threshold
        if steps == 0:
            residual = event.step - steps + self.mouse_wheel_step_residual
            adjust = copysign(floor(abs(residual)), residual)
            steps += adjust
            steps = round(steps)
            self.mouse_wheel_step_residual = residual - adjust  # update the residual

        # move slice by steps
        self.iv.move_slice(steps)

    def on_resize(self, event: mpl.backend_bases.ResizeEvent):
        """ Changes the dpi of the figure, so that the resolution does not change.
        Aims to increase performance on large window sizes.
        *** not done, not used ***

        @param event: the mpl.backend_bases.ResizeEvent to handle
        """
        size_inches = self.iv.canvas.figure.get_size_inches()
        new_dpi = self.iv.max_resolution // min(size_inches)
        # print(size_inches, new_dpi)
        # self.iv.canvas.figure.set_dpi(new_dpi)
        # self.iv.canvas.figure.set_size_inches(size_inches)
        # self.iv.canvas.draw()
        # self.iv.dpi = new_dpi
        # self.iv.redraw_slice()


if __name__ == '__main__':
    import sys

    if QtWidgets.QApplication.instance() is None:
        app = QtWidgets.QApplication(sys.argv)

    w = IMIImageViewer(None)
    w.show()
    sys.exit(QtWidgets.QApplication.instance().exec_())
