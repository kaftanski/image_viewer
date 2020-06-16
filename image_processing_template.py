import sys
import threading
import image_viewing as vis
import SimpleITK as sitk
from PyQt5.QtWidgets import QApplication

from lightweight_viewer import LightWeightViewer


class WindowManager:
    def __init__(self):
        if QApplication.instance():
            self.app = QApplication.instance()
        else:
            self.app = QApplication(sys.argv)

        self.windows = []
        #sys.exit(self.app.exec_())

    def show_image(self, image, window_title, blocking=True):
        exec_code = self._show_image(image, window_title, blocking=blocking)
        return exec_code

    def show_image_with_mask(self, image, mask, window_title, blocking=True):
        exec_code = self._show_image(image, window_title, blocking=blocking, mask=mask)
        return exec_code

    def show_and_return_markers(self, image, window_title):
        # window should be blocking by nature of the application
        markers = self._show_image(image, window_title, blocking=True, return_markers=True)
        return markers

    def _show_image(self, image, window_title, blocking=True, mask=None, return_markers=False):
        if not blocking and return_markers:
            raise ValueError('A viewer with non-blocking behaviour cannot return markers for further processing.')

        w = LightWeightViewer(image, window_title, mask=mask)
        w.show()
        self.windows.append(w)

        if mask is not None:
            w.image_viewer.add_mask(mask)
        # if blocking:
        #     threading.current_thread().join()
        if return_markers:
            return w.get_markers_for_region_growing()
        else:
            return 0

    def run(self):
        while True:
            if not vis.request_queue.empty():
                request = vis.request_queue.get()
                result = self._show_image(request.image, request.window_title, request.blocking, request.mask, request.return_markers)
                request.set_result(result)




# gets executed on worker thread
def process():
    # load example image from argv
    input_image = sitk.ReadImage(sys.argv[1])
    vis.show_image(input_image, 'Input Image (non-blocking window)', blocking=False)

    # pre-processing
    print('smoothing')
    smoothed = sitk.DiscreteGaussian(input_image, 2)
    vis.show_image(smoothed, 'Smoothed Image')

    # get seedpoints
    seeds = vis.show_and_return_markers(smoothed, 'Set Seedpoints')
    print(seeds)

    # segmentation
    binary_segmentation = sitk.ConnectedThreshold(smoothed, seeds, lower=650, upper=1000)

    # visualize result
    vis.show_image_with_mask(input_image, binary_segmentation, 'Result')


if __name__ == '__main__':
    window_manager = WindowManager()
    worker_thread = threading.Thread(target=process, name='worker')
    worker_thread.start()
    window_manager.run()
