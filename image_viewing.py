from PyQt5.QtWidgets import QApplication
import threading, queue

from lightweight_viewer import LightWeightViewer
import SimpleITK as sitk
import sys


_open_viewers = []

request_queue = queue.Queue(maxsize=5)


class ViewRequest:
    def __init__(self, image, window_title, blocking=True, mask=None, return_markers=False):
        self.image = image
        self.window_title = window_title
        self.blocking = blocking
        self.mask = mask
        self.return_markers = return_markers

        self.result = None

    def set_result(self, res):
        self.result = res

    def get_result(self):
        return self.result


def show_image(image, window_title, blocking=True):
    request = ViewRequest(image, window_title, blocking)
    if blocking:
        request_queue.join()
        result = request.get_result()
    else:
        result = 0

    return result
    # else:
    #     gui_thread = threading.Thread(target=_show_image, args=(image, window_title))
    #     gui_thread.start()
    #     return 0


def show_image_with_mask(image, mask, window_title, blocking=True):
    exec_code = _show_image(image, window_title, blocking=blocking, mask=mask)
    return exec_code


def show_and_return_markers(image, window_title):
    # window should be blocking by nature of the application
    markers = _show_image(image, window_title, blocking=True, return_markers=True)
    return markers


def _show_image(image, window_title, blocking=True, mask=None, return_markers=False):
    if not blocking and return_markers:
        raise ValueError('A viewer with non-blocking behaviour cannot return markers for further processing.')
    # print(threading.current_thread())
    # if QApplication.instance() is not None:
    #     app = QApplication.instance()
    # else:
    #     app = QApplication(sys.argv)
    app = QApplication(sys.argv)
    w = LightWeightViewer(image, window_title, mask=mask)
    global _open_viewers
    _open_viewers.append(w)
    w.show()

    # for viewer in _open_viewers:
    #     viewer.show()

    if blocking:
        exec_code = app.exec_()
    else:
        exec_code = 0

    if mask is not None:
        w.image_viewer.add_mask(mask)

    if return_markers:
        return w.get_markers_for_region_growing()
    else:
        return exec_code


if __name__ == '__main__':
    for m in show_and_return_markers(None, 'test'):
        print(m)
