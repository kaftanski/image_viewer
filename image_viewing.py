from PyQt5.QtWidgets import QApplication

from lightweight_viewer import LightWeightViewer, ImageMask
import SimpleITK as sitk
import sys


viewer_queue = []


def show_image(image, window_title, blocking=True):
    exec_code = _show_image(image, window_title, blocking=blocking)
    return exec_code


def show_image_with_mask(image, mask, window_title, color='b', blocking=True):
    im = ImageMask(mask, color=color)
    exec_code = _show_image(image, window_title, mask=im, blocking=blocking)
    return exec_code


def show_and_return_markers(image, window_title):
    # blocking should always be true in this application
    markers = _show_image(image, window_title, return_markers=True, blocking=True)
    return markers


def _show_image(image, window_title, blocking=True, mask=None, return_markers=False):
    # get the active Qt application or create a new one
    if QApplication.instance() is not None:
        app = QApplication.instance()
    else:
        app = QApplication(sys.argv)

    # construct the new viewer with the input image
    widget = LightWeightViewer(image, window_title, mask=mask)
    viewer_queue.append(widget)

    if blocking:
        # start app and show all previously invisible viewers
        for v in viewer_queue:
            v.show()
        exec_code = app.exec_()  # start the Qt event loop (blocking all further execution)

        # remove the viewers already shown
        viewer_queue.clear()
    else:
        # nothing wrong because no execution of event loop
        exec_code = 0

    if return_markers:
        # get the user input of the new viewer
        return widget.get_markers_for_region_growing()
    else:
        return exec_code


if __name__ == '__main__':
    # print(show_image_with_mask(None, sitk.ReadImage('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.OT_reg.nii.gz')))
    for m in show_and_return_markers(None, 'test'):
        print(m)
