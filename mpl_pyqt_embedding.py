from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton
import sys
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import itk

from mpl_interactive_view import show_in_interactive_figure


class Window(QMainWindow):
    def __init__(self):
        super().__init__()

        title = "MPL in PyQt5"
        top = 400
        left = 400
        width = 900
        height = 500

        self.setWindowTitle(title)
        self.setGeometry(top, left, width, height)
        self.my_ui()

    def my_ui(self, width=8, height=4):
        canvas = Canvas(self, 8, 4)
        canvas.move(0, 0)

        button = QPushButton('Clicky', parent=self)
        button.move(100, 450)

        button = QPushButton('Clacky', parent=self)
        button.move(250, 450)


class Canvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=5, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        fig = show_in_interactive_figure(itk.imread('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.MR_DWI_reg.nii.gz'), fig)
        #self.plot()

    def plot(self):
        x = np.arange(10)
        y = x**2
        self.axes.plot(x, y)


app = QApplication(sys.argv)
window = Window()
window.show()
app.exec()
