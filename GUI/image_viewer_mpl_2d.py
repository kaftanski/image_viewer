from typing import Union, Sequence, Tuple

import matplotlib as mpl
import numpy as np
from matplotlib import image as mpimg
from matplotlib import pyplot as plt
from matplotlib.backend_bases import MouseButton, cursors


# mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['image.cmap'] = 'gray'


class Image:
    def __init__(self, image_array: np.ndarray, spacing: Sequence[float] = None, origin: Sequence[Union[int, float]] = None):
        # ensure 2D or 3D image
        image_dim = len(image_array.shape)
        assert 2 <= image_dim <= 3

        self.image_array = image_array
        self.spacing = spacing if spacing is not None else np.ones(image_dim)
        self.origin = origin if origin is not None else np.zeros(image_dim)

    def GetImageArray(self):
        """
        @return: the numpy array containing the image intensities
        """
        return self.image_array

    def GetSpacing(self) -> Sequence[float]:
        """
        @return: the image spacing
        """
        return self.spacing

    def GetOrigin(self) -> Sequence[Union[int, float]]:
        """
        @return: the image origin
        """
        return self.origin

    def GetSize(self) -> Tuple[int]:
        """
        @return: the dimensions of the image
        """
        return self.image_array.shape

    def GetPixel(self, x: int, y: int, z: int = None) -> Union[int, float]:
        if z is None:
            assert len(self.GetSize()) == 2, 'Error: Trying to access a 3D image with only 2 indices. ' \
                                             'If needed, use numpy-like indexing with []-brackets.'
            return self.image_array[x, y]
        else:
            assert len(self.GetSize()) == 3, 'Error: Trying to access a 2D image with 3 indices.'
            return self.image_array[x, y, z]

    def SetPixel(self, value: Union[int, float], x: int, y: int, z: int = None):
        if z is None:
            assert len(self.GetSize()) == 2, 'Error: Trying to access a 3D image with only 2 indices. ' \
                                             'If needed, use numpy-like indexing with []-brackets.'
            self.image_array[x, y] = value
        else:
            assert len(self.GetSize()) == 3, 'Error: Trying to access a 2D image with 3 indices.'
            self.image_array[x, y, z] = value

    def __getitem__(self, item) -> np.ndarray:
        return self.image_array[item]

    def __setitem__(self, key, value):
        self.image_array[key] = value


class MPLImageViewer:
    def __init__(self, image_data: Image = None):
        if image_data is None:
            # create an empty 3d image
            image_data = Image(np.zeros((1, 1, 1)))

        # init attributes
        self.image = image_data  # the SimpleITK image in the viewer
        self.image_array = image_data.GetImageArray()  # the image as np.array
        self.im_plot = None

        self.current_window = 0
        self.current_level = 0

        self.greyval_range = [0, 0]

        self.masks = {}  # dict with mask as key and its plot as value

        initial_canvas_size_x = 5  # size in pixels. corresponds to a size of 5x5 inches with 100 dpi
        initial_canvas_size_y = 5
        self.max_resolution = 500
        self.dpi = self.max_resolution // initial_canvas_size_x

        self.canvas = plt.figure(figsize=(initial_canvas_size_x, initial_canvas_size_y), dpi=self.dpi, facecolor='black').canvas

        self.ax = self.canvas.figure.subplots()
        self.ax.axis('off')
        self.canvas.figure.subplots_adjust(0, 0, 1, 1)

        # interaction
        self.toolbar = self.canvas.toolbar  # use the built-in toolbar, which is added on the plt.figure() call
        self.toolbar.hide()

        self.window_level_start_position = None
        self.init_interaction()

        # show the viewer
        self.init_image()
        plt.show()

    def init_image(self):
        """ Initial call after the image has been set or reset.

        Constructs the image array and plots the middle slice of the image with window/level spanning the whole
        grey value range of the image.
        """
        # DON'T swap axes, in order to keep the fastest dimension of the np array the first
        # -> z-slices are the most performant
        self.image_array = self.image.GetImageArray()

        # redraw if there is already a plot
        if self.im_plot is not None:
            self.im_plot.remove()

        aspect_ratio = (self.image.GetSpacing()[1] * self.image.GetSize()[1]) / (self.image.GetSpacing()[0] * self.image.GetSize()[0])
        self.im_plot = self.ax.imshow(self.image_array, aspect=aspect_ratio)

        # set the window/level to span the whole greyvalue range of the image
        min_grey = self.image_array.min()
        max_grey = self.image_array.max()
        self.greyval_range = [min_grey, max_grey]
        self.reset_window_level()

    def init_interaction(self):
        self.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.canvas.mpl_connect('button_press_event', self.handle_mouse_button_down)
        self.canvas.mpl_connect('button_release_event', self.handle_mouse_button_up)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_motion)

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
        self.canvas.draw()
        print('Window: {}, Level: {}'.format(self.current_window, self.current_level))

    def reset_window_level(self):
        """ Set the window/level to span the entire grey value range of the image.
        """
        # reset to cover the entire image range
        window = self.greyval_range[1] - self.greyval_range[0]
        level = self.greyval_range[0] + window / 2
        self.set_window_level(window, level)

    def on_key_press(self, event):
        if event.key == 'r':
            self.reset_window_level()
            
    def handle_mouse_button_down(self, event: mpl.backend_bases.MouseEvent):
        """ Handles mouse button down events by distributing event to
        self.on_left_button_down or self.on_right_button_down respectively.

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
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
        pass

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
        self.toolbar.pan()

        # update cursor manually (otherwise this is only done after the mouse is first moved)
        if self.toolbar.mode == "pan/zoom":
            self.toolbar.set_cursor(cursors.MOVE)  # grabbing hand
        else:
            self.toolbar.set_cursor(cursors.POINTER)  # normal cursor

    def on_right_button_down(self, event: mpl.backend_bases.MouseEvent):
        """ Left mouse button down event handler:
             starts the window levelling.
             Moving the mouse while the right mouse button is pressed does window levelling.
             Movement in x-direction modifies the window, y-direction for level.

        @param event: the mpl.backend_bases.MouseEvent to handle
        """
        # if panning / zooming is not activated
        if self.toolbar.mode == '':
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
            new_window = self.current_window + dx * gain
            new_level = self.current_level + dy * gain/2

            self.set_window_level(new_window, new_level)


if __name__ == '__main__':
    MPLImageViewer(Image(mpimg.imread('/home/paul/Documents/imi_projects/MBV/MIPImages/lena.png')[:, :, 0]*255))
