from PyQt5.QtWidgets import QApplication

from lightweight_viewer import LightWeightViewer
import SimpleITK as sitk
import sys


def show_image(image, window_title, blocking=True):
    exec_code = _show_image(image, window_title)
    return exec_code


def show_image_with_mask(image, mask, window_title):
    exec_code = _show_image(image, window_title, mask=mask)
    return exec_code


def show_and_return_markers(image, window_title):
    markers = _show_image(image, window_title, return_markers=True)
    return list(list(int(coord) for coord in m.world_position) for m in markers)  # TODO: markers have to be int for region growing


def _show_image(image, window_title, blocking=True, mask=None, return_markers=False):
    markers = []

    app = QApplication(sys.argv)
    w = LightWeightViewer(image, window_title, mask=mask, p_markers=markers)
    exec_code = app.exec_()
    # TODO: non blocking viewers!
    if mask is not None:
        w.image_viewer.add_mask(mask)

    if return_markers:
        return markers
    else:
        return exec_code


if __name__ == '__main__':
    # print(show_image_with_mask(None, sitk.ReadImage('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.OT_reg.nii.gz')))
    for m in show_and_return_markers(None, 'test'):
        print(m)
