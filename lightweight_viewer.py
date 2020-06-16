from math import floor, copysign
from time import time

import SimpleITK as sitk
import numpy as np
from PyQt5.QtWidgets import QAction
import os

import matplotlib as mpl
from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.backend_bases import MouseButton

mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['image.cmap'] = 'gray'


SLICE_ORIENTATION_XY = 2
SLICE_ORIENTATION_XZ = 1
SLICE_ORIENTATION_YZ = 0

SLICE_ORIENTATION = {'xy': 2, 'xz': 1, 'yz': 0}

COLORS = {
    'red': [1, 0, 0],
    'green': [0, 1, 0],
    'blue': [0, 0, 1]
}


PARAMS = {
    'marker_radius': 3
}


class LightWeightViewer(QtWidgets.QMainWindow):
    def __init__(self, image, title='IMI Image Viewer', mask=None):
        super(LightWeightViewer, self).__init__()

        self.image_viewer = ImageViewer(image, parent=self)
        if mask is not None:
            self.image_viewer.add_mask(mask)

        self.setCentralWidget(self.image_viewer)
        # self.addToolBar(NavigationToolbar(self.image_viewer.canvas, self))
        self.init_menu_bar()
        self.setWindowTitle(title)
        self.show()

    def get_markers_for_region_growing(self):
        return list(list(int(coord) for coord in m.pixel_position) for m in self.image_viewer.markers)

    def init_menu_bar(self):
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

        clear_markers_action = QAction('&Clear all markers', self)
        clear_markers_action.setStatusTip('Clear all markers from all images')
        clear_markers_action.triggered.connect(lambda: self.image_viewer.remove_markers(only_last=False))

        remove_last_marker_action = QAction('Remove &last marker', self)
        remove_last_marker_action.setStatusTip('Remove the last placed marker')
        remove_last_marker_action.setShortcut('Ctrl+Z')
        remove_last_marker_action.triggered.connect(lambda: self.image_viewer.remove_markers(only_last=True))

        marker_submenu = view_menu.addMenu('&Marker')
        marker_submenu.addActions([clear_markers_action, remove_last_marker_action])

        zoom_action = QAction('&Zoom...', self)
        zoom_action.setStatusTip('Choose a zoom factor')
        zoom_action.triggered.connect(self.zoom_dialog)

        reset_wl_action = QAction('&Reset Window / Level', self)
        reset_wl_action.setStatusTip('Reset window-levelling to span the entire image range')
        reset_wl_action.triggered.connect(self.image_viewer.reset_window_level)

        view_menu.addActions([zoom_action, reset_wl_action])

    # def keyPressEvent(self, e):
    #     if e.text() == 'x':
    #         self.image_viewer.change_orientation(SLICE_ORIENTATION['yz'])
    #     elif e.text() == 'y':
    #         self.image_viewer.change_orientation(SLICE_ORIENTATION['xz'])
    #     elif e.text() == 'z':
    #         self.image_viewer.change_orientation(SLICE_ORIENTATION['xy'])

    def zoom_dialog(self):
        new_zoom_factor, ok = QtWidgets.QInputDialog.getInt(self, 'Zoom', 'Percentage:', min=100, step=20,
                                                        value=self.image_viewer.current_zoom)
        if ok:
            self.image_viewer.zoom(percent=new_zoom_factor)

    def file_dialog(self, caption, filter, save=False):
        if save:
            file_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, caption=caption,
                                                                 directory=os.getcwd(),
                                                                 filter=filter)
        else:
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, caption=caption, directory=os.getcwd(), filter=filter)
        return file_path

    def load_image(self):
        file_path = self.file_dialog(caption='Load Image', filter='Image (*.nii.gz)', save=False)

        if file_path != '':
            print('Loading image ' + file_path)

            image = sitk.ReadImage(file_path)
            self.image_viewer.set_image(image)

    def load_mask(self):
        file_path = self.file_dialog(caption='Load Image Mask', filter='Image (*.nii.gz)', save=False)

        if file_path != '':
            print('Loading mask ' + file_path)

            image = sitk.ReadImage(file_path)
            mask = ImageMask(image, color='r')
            self.image_viewer.add_mask(mask)

    def save_current_view(self):
        file_path = self.file_dialog('Save File', 'Image (*.png, *.PNG)', save=True)

        if file_path != '':
            if file_path[-4:].lower() != '.png':
                file_path += '.png'

            print('Saving view to {}'.format(file_path))
            self.image_viewer.canvas.print_figure(filename=file_path, dpi=200, bbox_inches=0)


class ImageViewer(QtWidgets.QWidget):
    def __init__(self, image_data: sitk.Image = None, parent=None):
        super(ImageViewer, self).__init__(parent)
        self.viewportsize_x = 800
        self.viewportsize_y = 800
        self.masks = {}  # dict with mask as key and its plot as value
        self.markers = []
        self.orientation = SLICE_ORIENTATION_XY

        if image_data is None:
            # load example image
            image_data = sitk.ReadImage('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.MR_DWI_reg.nii.gz')

        self.image = image_data
        self.image_array = None
        self.current_slice = 0

        self.im_plot = None

        self.current_window = 0
        self.current_level = 0

        self.greyval_range = [0, 0]
        self.current_zoom = 100

        # init the MPL canvas
        self.canvas = FigureCanvas(Figure(figsize=(5, 5), facecolor='white'))
        self.canvas.setParent(self)
        self.ax = self.canvas.figure.subplots()
        self.ax.axis('off')
        self.canvas.figure.subplots_adjust(0, 0, 1, 1)

        self.marker_plot = self.ax.scatter([], [], c=ImageMarker.STANDARD_COLOR)

        # define a slider for image slicing
        self.slice_slider = QtWidgets.QSlider(orientation=QtCore.Qt.Vertical, parent=self)
        self.slice_slider.valueChanged.connect(self.update_slice)

        # define a text field to display image coordinates
        bottom_bar_layout = QtWidgets.QHBoxLayout()
        self.zoom_label = QtWidgets.QLabel('{} %'.format(self.current_zoom))
        self.pixel_info_label = PixelInfoQLabel(parent=self)
        bottom_bar_layout.addWidget(self.zoom_label)
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

        self.init_image()

        self.interactor_style = ImageViewerInteractor(self)

        self.show()

        # self.add_mask(sitk.ReadImage('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.OT_reg.nii.gz'))
        # self.add_mask(sitk.ReadImage('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/02/VSD.Brain.02.O.OT_reg.nii.gz'), color='g')

    def init_image(self):
        # DON'T swap axes, in order to keep the fastest dimension of the np array the first
        # -> z-slices are the most performant
        self.image_array = sitk.GetArrayFromImage(self.image)

        self.current_slice = self.get_slice_range()[1] // 2

        if self.im_plot is not None:
            self.im_plot.remove()

        self.im_plot = self.ax.imshow(self.image_array[index_compatibility(get_3d_plane_index(self.current_slice, self.orientation))],
                                      aspect=get_aspect_ratio_for_plane(self.image.GetSpacing(), self.orientation, self.image.GetSize()))

        min_grey = self.image_array.min()
        max_grey = self.image_array.max()
        self.greyval_range = [min_grey, max_grey]
        self.reset_window_level()

        initial_coords = [0, 0]
        initial_coords.insert(self.orientation, self.current_slice)
        self.show_pixel_info(initial_coords)

        self.update_slider()

    def change_orientation(self, orientation):
        if orientation not in SLICE_ORIENTATION.keys() and orientation not in SLICE_ORIENTATION.values():
            print('Cannot change viewer orientation to {}. Use one of {} or {} respectively.'.
                  format(orientation, SLICE_ORIENTATION.keys(), SLICE_ORIENTATION.values()))

        if isinstance(orientation, str):
            orientation = SLICE_ORIENTATION[orientation]

        if orientation == self.orientation:
            return

        self.orientation = orientation
        self.update_slice(self.current_slice)
        self.update_slider()

    def set_window_level(self, window, level):
        if window < 0:
            # set a minimum level (could be even smaller in case of normalized images with values in [0,1])
            window = 1e-5

        self.current_window = window
        self.current_level = level

        half_window = self.current_window / 2
        self.im_plot.set_clim(vmin=self.current_level - half_window, vmax=self.current_level + half_window)
        # self.ax.draw_artist(self.im_plot)
        self.canvas.draw()  # TODO: only redraw the image, not markers or masks
        print('Window: {}, Level: {}'.format(self.current_window, self.current_level))

    def reset_window_level(self):
        # reset to cover the entire image range
        window = self.greyval_range[1] - self.greyval_range[0]
        level = self.greyval_range[0] + window / 2
        self.set_window_level(window, level)

    def update_slice(self, new_slice):
        # todo: Blitting
        lower, upper = self.get_slice_range()
        if not lower <= new_slice < upper:
            new_slice = lower if new_slice < lower else upper - 1

        self.current_slice = new_slice

        img_slice = self.image_array[index_compatibility(get_3d_plane_index(slice_index=new_slice, orientation=self.orientation))]

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

        self.pixel_info_label.set_coordinate(self.orientation, new_slice)
        self.show_pixel_info(self.pixel_info_label.coords)  # TODO: move this to label class?

        # reset zoom factor
        self.zoom_label.setText('100 %')
        self.current_zoom = 100

        self.update_masks()
        self.scatter_markers()

        self.canvas.draw()

    def zoom(self, percent, center_x=None, center_y=None):
        if percent < 100:
            return

        self.current_zoom = percent

        plane_dimensions = list(self.image.GetSize()[d] for d in range(3) if d != self.orientation)

        if center_x is None:
            center_x = (plane_dimensions[0] - 1) / 2
        if center_y is None:
            center_y = (plane_dimensions[1] - 1) / 2

        half_range_x = plane_dimensions[0] / 2 / (percent / 100)
        half_range_y = plane_dimensions[1] / 2 / (percent / 100)

        x_limits = [center_x - half_range_x, center_x + half_range_x]
        y_limits = [center_y - half_range_y, center_y + half_range_y]

        # adjust limits to keep the image in the figure
        for dim, limits in enumerate([x_limits, y_limits]):
            shift = 0
            if limits[0] < 0:
                shift = -limits[0]
            elif limits[1] > plane_dimensions[dim] - 1:
                shift = ((plane_dimensions[dim] - 1) - limits[1])

            limits[0] += shift
            limits[1] += shift

        self.ax.set_xlim(*x_limits)
        self.ax.set_ylim(*y_limits)

        self.zoom_label.setText('{} %'.format(percent))

        self.canvas.draw()

    def get_slice_range(self):
        return 0, self.image.GetSize()[self.orientation]

    def update_slider(self):
        self.slice_slider.setRange(self.get_slice_range()[0], self.get_slice_range()[1]-1)
        self.slice_slider.setValue(self.current_slice)

    def move_slice(self, delta):
        # moves the slider, which in turn changes the slice
        self.slice_slider.setValue(self.current_slice + delta)  # TODO: kind of a workaround

    def add_marker(self, position):
        assert len(position) == 3
        marker = ImageMarker(position)
        self.markers.append(marker)
        self.scatter_markers()
        self.canvas.draw()

    def remove_markers(self, only_last=False):
        if len(self.markers) == 0:
            return

        if only_last:
            self.markers.pop()
        else:
            self.markers.clear()

        self.scatter_markers()
        self.canvas.draw()

    def scatter_markers(self):
        markers_in_slice = []
        for m in self.markers:
            pixel_pos = np.array(m.pixel_position)
            if abs(pixel_pos[self.orientation] - self.current_slice) < 0.5:
                markers_in_slice.append(pixel_pos[np.arange(3) != self.orientation])

        if len(markers_in_slice) != 0:
            np_markers = np.stack(markers_in_slice, axis=1)
            self.marker_plot.set_offsets(np_markers.swapaxes(0, 1))
        else:
            self.marker_plot.remove()
            self.marker_plot = self.ax.scatter([], [], c=ImageMarker.STANDARD_COLOR)

    def add_mask(self, mask):
        if not compatible_metadata(self.image, mask.itk_mask):
            return

        self.masks[mask] = None
        self.update_masks()
        self.canvas.draw()

    def update_masks(self):
        self.clear_masks()
        for m in self.masks.keys():
            self.masks[m] = add_mask_to_image(self.ax, m.get_slice(self.current_slice, self.orientation),
                                              aspect=get_aspect_ratio_for_plane(m.get_spacing(), self.orientation, self.image.GetSize()),
                                              alpha=m.alpha, color=m.color)

    def set_image(self, image):
        self.image = image

        # remove old masks
        self.clear_masks()
        self.masks.clear()

        # remove old markers
        self.markers.clear()
        self.scatter_markers()

        # draw image and re-init everything else
        self.init_image()

    def clear_masks(self):
        for m in self.masks.keys():
            try:
                # remove the existing plots
                self.masks[m].remove()
            except AttributeError:
                # if the mask was not there in the previous slice it is None here
                pass

    def show_pixel_info(self, pixel_coords):
        if pixel_coords is None:
            self.pixel_info_label.clear()
        else:
            x, y, z = [int(ind) for ind in pixel_coords]  # attention to different indexing of np array from itk image
            try:
                self.pixel_info_label.set_values(x, y, z, self.image.GetPixel(x, y, z))
            except (RuntimeError, TypeError):
                # gets thrown if x, y, z are out of bounds of the image (this happens on the edges of the figure)
                return


class PixelInfoQLabel(QtWidgets.QLabel):
    def __init__(self, parent=None):
        super(PixelInfoQLabel, self).__init__(text='', parent=parent)

        self.coords = [0, 0, 0]
        self.intensity = 0

    def set_values(self, x, y, z, intensity):
        self.coords = [x, y, z]
        self.intensity = intensity

        self.update_text()

    def set_coordinate(self, dim, coord):
        self.coords[dim] = coord
        self.update_text()

    def set_intensity(self, intensity):
        self.intensity = intensity
        self.update_text()

    def update_text(self):
        self.setText('({}, {}, {}): {}'.format(*self.coords, self.intensity))


class ImageViewerInteractor:
    def __init__(self, image_viewer: ImageViewer):
        self.iv = image_viewer

        self.iv.canvas.mpl_connect('button_press_event', self.handle_mouse_button_down)
        self.iv.canvas.mpl_connect('button_release_event', self.handle_mouse_button_up)
        self.iv.canvas.mpl_connect('motion_notify_event', self.on_mouse_motion)
        self.iv.canvas.mpl_connect('scroll_event', self.on_mousewheel_event)
        self.iv.canvas.mpl_connect('key_press_event', self.on_key_press)

        self.last_mouse_position_in_figure = [0, 0]
        self.window_level_start_position = None
        self.mouse_wheel_step_residual = 0  # accumulating mousewheel movements, if they dont add up to a whole step

    def on_key_press(self, event):
        print(event.key)
        if event.key == 'x':
            self.iv.change_orientation(SLICE_ORIENTATION['yz'])
        elif event.key == 'y':
            self.iv.change_orientation(SLICE_ORIENTATION['xz'])
        elif event.key == 'z':
            self.iv.change_orientation(SLICE_ORIENTATION['xy'])

    def handle_mouse_button_down(self, event):
        print(event.button)
        if event.button == MouseButton.LEFT:
            self.on_left_button_down(event)
        elif event.button == MouseButton.RIGHT:
            self.on_right_button_down(event)

    def handle_mouse_button_up(self, event):
        print(event.button, 'up')
        if event.button == MouseButton.LEFT:
            self.on_left_button_up(event)
        elif event.button == MouseButton.RIGHT:
            self.on_right_button_up(event)

    def on_mouse_motion(self, event):
        if event.xdata is not None and event.ydata is not None:
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

    def on_mousewheel_event(self, event):
        steps = round(event.step)

        # in case of high resolution touchpads/mice: accumulate values under the step threshold
        if steps == 0:
            residual = event.step - steps + self.mouse_wheel_step_residual
            adjust = copysign(floor(abs(residual)), residual)
            steps += adjust
            self.mouse_wheel_step_residual = residual - adjust  # update the residual

        if event.key == 'control':
            # if ctrl key is pressed while scrolling: zoom image

            if event.xdata is None or event.ydata is None:
                # choose the figure center as zoom center if mouse outside of figure
                center_x = self.last_mouse_position_in_figure[0]
                center_y = self.last_mouse_position_in_figure[1]
            else:
                center_x = event.xdata
                center_y = event.ydata

            # zooming steps * 20 %
            new_zoom = self.iv.current_zoom + steps * 20
            self.iv.zoom(percent=new_zoom, center_x=center_x, center_y=center_y)
        else:
            # move slice by steps
            self.iv.move_slice(steps)
            self.iv.pixel_info_label.set_coordinate(self.iv.orientation, self.iv.current_slice)
            self.iv.show_pixel_info(self.iv.pixel_info_label.coords)

    def on_left_button_down(self, event):
        if event.xdata is not None and event.ydata is not None:
            # click inside of figure
            event_data_position = [event.xdata, event.ydata]
            event_data_position.insert(self.iv.orientation, self.iv.current_slice)
            print('adding marker at {}'.format(event_data_position))
            self.iv.add_marker(event_data_position)

    def on_left_button_up(self, event):
        pass

    def on_right_button_down(self, event):
        # start window levelling
        self.window_level_start_position = [event.x, event.y]

    def on_right_button_up(self, event):
        self.window_level_start_position = None


class ImageMask:
    def __init__(self, itk_binary_image, alpha=0.3, color='r'):
        self.itk_mask = itk_binary_image
        self.mask_as_array = sitk.GetArrayFromImage(self.itk_mask)
        self.alpha = alpha
        self.color = color

    def get_slice(self, slice_index, orientation):
        return self.mask_as_array[index_compatibility(get_3d_plane_index(slice_index, orientation))]

    def get_spacing(self):
        return self.itk_mask.GetSpacing()


class ImageMarker:
    STANDARD_COLOR = 'b'

    def __init__(self, pixel_position):
        assert len(pixel_position) == 3
        self.pixel_position = pixel_position

    def __str__(self):
        return 'ImageMarker at pixel position ({}, {}, {})'.format(*self.pixel_position)


def get_3d_plane_index(slice_index, orientation):
    index = [slice(0, None), slice(0, None)]  # get the whole dimensions of planar axes
    index.insert(orientation, slice_index)
    return tuple(index)


def get_aspect_ratio_for_plane(spacing, orientation, image_dimensions):
    dims = list(d for d in range(3) if d != orientation)
    ratio = spacing[dims[1]] * image_dimensions[1] / (spacing[dims[0]] * image_dimensions[0])
    return ratio


def compatible_metadata(image1: sitk.Image, image2: sitk.Image, check_size=True, check_spacing=True, check_origin=True):
    all_parameters_equal = True
    tolerance = 1e-4

    if check_size:
        size1 = image1.GetSize()
        size2 = image2.GetSize()
        if size1 != size2:
            all_parameters_equal = False
            print(f'Images do not have the same size ({size1} != {size2})')

    if check_spacing:
        spacing1 = image1.GetSpacing()
        spacing2 = image2.GetSpacing()
        if any(list(abs(s1-s2) > tolerance for s1, s2 in zip(spacing1, spacing2))):
            all_parameters_equal = False
            print(f'Images do not have the same spacing ({spacing1} != {spacing2})')

    if check_origin:
        origin1 = image1.GetOrigin()
        origin2 = image2.GetOrigin()
        if any(list(abs(o1-o2) > tolerance for o1, o2 in zip(origin1, origin2))):
            all_parameters_equal = False
            print(f'Images do not have the same origin ({origin1} != {origin2})')

    return all_parameters_equal


def index_compatibility(index):
    # change to the other method (sitk uses x, y, z, numpy uses z, y, x
    assert len(index) == 3
    return index[-1::-1]


def pix2world(index, spacing):
    assert len(index) == len(spacing)
    return tuple(i*s for i, s in zip(index, spacing))


def world2pix(index, spacing):
    assert len(index) == len(spacing)
    return tuple(i/s for i, s in zip(index, spacing))


def add_mask_to_image(ax, mask, aspect, alpha=0.3, color='r'):
    """ can be multilabel, if mask is given as (x, y, l) with l different binary images """
    if mask.sum() == 0:
        return None

    cmap = ListedColormap([color])
    mask = np.ma.masked_where(mask == 0, mask)
    plot = ax.matshow(mask, aspect=aspect, cmap=cmap, alpha=alpha)

    return plot


if __name__ == '__main__':
    import sys

    if QtWidgets.QApplication.instance() is None:
        app = QtWidgets.QApplication(sys.argv)
    # app = QtWidgets.QApplication(sys.argv)
    w = LightWeightViewer(None)
    w.show()
    # print('test')
    # for i in range(1000000):
    #     print(i)
    sys.exit(QtWidgets.QApplication.instance().exec_())
