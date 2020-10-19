import atexit
import sys
from typing import List, Union

import SimpleITK as sitk
from PyQt5.QtWidgets import QApplication

from viewer_utils import ImageMask
from imi_image_viewer import LightWeightViewer

_viewer_queue = []


def show_image(image: sitk.Image, window_title: str = '', blocking: bool = True) -> int:
    """ Display an image.

    @param image: the image to show
    @param window_title: optional window title of the viewer application
    @param blocking: if set to true, the viewer will be shown and following code is continued after the viewer is closed
    @return: execution code of the viewer application. this equals 0 if the application terminated without error.
    """
    exec_code = _show_image(image, window_title, blocking=blocking)
    return exec_code


def show_image_with_mask(image: sitk.Image, mask_image: sitk.Image, window_title: str = '', color: str = 'b', blocking: bool=True) -> int:
    """ Display an image with a mask overlay

    @param image: the image to show
    @param mask_image: the mask image to show on top of the image (should have the same size, spacing and origin as the image)
    @param window_title: optional window title of the viewer application
    @param color: optionally specify the color of the mask (look up 'specifying colors in matplotlib' for different methods of defining a color)
    @param blocking: if set to true, the viewer will be shown and following code is continued after the viewer is closed
    @return: execution code of the viewer application. this equals 0 if the application terminated without error.
    """
    im = ImageMask(mask_image, color=color)
    exec_code = _show_image(image, window_title, mask=im, blocking=blocking)
    return exec_code


def show_and_return_markers(image: sitk.Image, window_title: str = '') -> List[List[int]]:
    """ Display an image and return the markers that have been manually placed in the viewer.

    @param image: the image to show
    @param window_title: optional window title of the viewer application
    @return: a list of markers (which are each a list containing 3 indices)
    """
    # blocking is always true in this use case
    markers = _show_image(image, window_title, return_markers=True, blocking=True)
    return markers


def _show_image(image: sitk.Image, window_title: str = '', blocking: bool = True, mask: ImageMask = None, return_markers: bool = False) -> Union[int, List[List[int]]]:
    """ Function to access all image viewing functionality (do not call this directly!)

    @param image: the image to show
    @param window_title: optional window title of the viewer application
    @param blocking: if set to true, the viewer will be shown and following code is continued after the viewer is closed
    @param mask: an optional image mask to show on top of the image
    @param return_markers: if true, a list of markers manually placed in the viewer, otherwise the execution code is returned.
    @return: execution code of the viewer application if return_markers is false, otherwise the markers are returned.
    """
    # get the active Qt application or create a new one
    if QApplication.instance() is None:
        app = QApplication(sys.argv)

    # construct the new viewer with the input image
    widget = LightWeightViewer(image, window_title, mask=mask)
    _viewer_queue.append(widget)

    if blocking:
        # start app and show all previously invisible viewers
        exec_code = _show_viewers_in_queue()
    else:
        # nothing wrong because no execution of event loop
        exec_code = 0

    if return_markers:
        # get the user input of the new viewer
        return widget.get_markers_for_region_growing()
    else:
        return exec_code


def _show_viewers_in_queue() -> int:
    """ Shows all viewers in _viewer_queue.

    @return: execution code of the viewer application. this equals 0 if the application terminated without error.
    """
    app = QApplication.instance()

    # show viewers, if app exists and there are viewers left to show
    if app is not None and len(_viewer_queue) > 0:
        for v in _viewer_queue:
            v.show()
        exec_code = app.exec_()  # start the Qt event loop (blocking all further execution)
        # remove the viewers already shown
        _viewer_queue.clear()
        return exec_code
    else:
        return 0


# register exit hook that shows all previously not shown viewers (happens when the last opened viewer is non blocking)
atexit.register(_show_viewers_in_queue)


if __name__ == '__main__':
    # print(show_image_with_mask(None, sitk.ReadImage('/home/paul/Documents/imi_projects/MBV/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.OT_reg.nii.gz')))
    for m in show_and_return_markers(None, 'test'):
        print(m)
