import itk
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib import colors as mpc

import numpy as np
mpl.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib.backend_bases import MouseButton


MIN_WINDOW = 10
MIN_LEVEL = -1000
MAX_LEVEL = 3000


class SliceHandler(object):
    def __init__(self, ax, img_array):
        self.ax = ax
        self.img_array = img_array
        _, _, self.slices = img_array.shape
        self.ind = self.slices // 2

        self.im = ax.imshow(self.img_array[:, :, self.ind], cmap='gray')

        image_limits = self.im.get_clim()
        self.window = image_limits[1] - image_limits[0]
        self.level = image_limits[0] + self.window / 2

        self.last_left_mouseclick = None

        self.update()
        self.update_window_level()

    def onscroll(self, event):
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def onmousedown(self, event):
        if event.button == MouseButton.LEFT:
            #print('lm click at ({}, {})'.format(event.x, event.y))
            self.last_left_mouseclick = [event.x, event.y]

    def onmouserelease(self, event):
        if event.button == MouseButton.LEFT:
            self.last_left_mouseclick = None

    def onmotion(self, event):
        #print('movement at ({}, {})'.format(event.xdata, event.ydata))
        #print(event.x, event.y)
        if self.last_left_mouseclick is None:
            return

        # move window by x movement
        # move level by y movement
        dx = event.x - self.last_left_mouseclick[0]
        dy = event.y - self.last_left_mouseclick[1]

        if event.x is not None:
            self.last_left_mouseclick[0] = event.x
        if event.y is not None:
            self.last_left_mouseclick[1] = event.y

        self.window += dx
        #  print(dy)
        self.level += dy

        if self.window < MIN_WINDOW:
            self.window = MIN_WINDOW

        if self.level < MIN_LEVEL:
            self.level = MIN_LEVEL
        elif self.level > MAX_LEVEL:
            self.level = MAX_LEVEL

        self.update_window_level()

    def update_window_level(self):
        half_window = self.window / 2
        self.im.set_clim(vmin=self.level - half_window, vmax=self.level + half_window)
        self.im.axes.figure.canvas.draw()
        print('Window: {}, Level: {}'.format(self.window, self.level))

    def update(self):
        self.im.set_data(self.img_array[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()


class Interactive3DFigure(Figure):
    def __init__(self, image_data):
        super(Interactive3DFigure, self).__init__(tight_layout=True)


def window_level(image, w, l):
    half_window = w / 2
    image.set_clim(vmin=l - half_window, vmax=l + half_window)
    image.axes.figure.canvas.draw()


def add_mask_image(ax, mask):
    colors = ('b', 'y', 'g', 'r', 'c')
    labels = np.unique(mask)

    for label in labels:
        if label == 0:
            continue

        cmap = mpc.ListedColormap([colors[label]])
        mask = np.ma.masked_where(mask == 0, mask)
        ax.matshow(mask, cmap=cmap, alpha=0.5)

    return ax


def show_in_interactive_figure(itk_image, fig):
    # fig = Figure(figsize=(5, 5), dpi=100)
    # ax = fig.add_subplot(111)
    ax = fig.get_axes()[0]
    handler = SliceHandler(ax, itk.array_view_from_image(itk_image))
    fig.canvas.mpl_connect('scroll_event', handler.onscroll)
    fig.canvas.mpl_connect('button_press_event', handler.onmousedown)
    fig.canvas.mpl_connect('button_release_event', handler.onmouserelease)
    fig.canvas.mpl_connect('motion_notify_event', handler.onmotion)
    ax.plot(105, 60, 'ro')
    return fig


if __name__ == '__main__':
    fig = plt.figure(figsize=(5, 5), dpi=100)
    ax = fig.add_subplot(111)
    test_image = itk.imread('/home/paul/Documents/imi_projects/MBV/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.MR_DWI_reg.nii.gz')
    fig = show_in_interactive_figure(test_image, fig)

    maskimage = itk.imread('/home/paul/Documents/imi_projects/MBV/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.OT_reg.nii.gz')
    maskimage = itk.array_view_from_image(maskimage)[60]
    # add_mask_image(ax, maskimage)
    plt.show()
