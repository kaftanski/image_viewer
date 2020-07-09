import SimpleITK as sitk
from PyQt5 import QtWidgets
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QAction

from vtk_viewer.vtk_image_viewer import ImageViewer


class ImageContainer:
    def __init__(self):
        self.image_map = {}

    def add_image(self, key, image):
        self.image_map[key] = image

    def get_image(self, key):
        return self.image_map[key]


class MainWindowHandler:
    def __init__(self, main_window):
        self.main_window = main_window
        self.image_container = ImageContainer()

        self.setup_connections()

    def setup_connections(self):
        pass

    def load_image(self):
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self.main_window, caption='Load Image', directory='/home/paul/Documents/imi_projects/MBV/MIPImages', filter='Image (*.nii.gz)')
        print('Loading image ' + file_path)
        try:
            file_name = file_path.split('/')[-1]
            image = sitk.ReadImage(file_path)
            self.image_container.add_image(file_name, image)
            self.main_window.add_image(image, file_name)
        except RuntimeError:
            print('No file selected')


class ProjectMainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.handler = MainWindowHandler(self)
        self.content = ProjectMainWindowContent(self)

        # central vtk frame
        self.setCentralWidget(self.content)

        # controls
        self.init_menuBar()
        self.statusBar()  # TODO: show xyz coordinates and greyvalue in real time here

        self.setGeometry(100, 100, 1200, 800)
        self.setWindowTitle('MIP Project UI')
        self.show()

        self.add_image(sitk.ReadImage(
            '/home/paul/Documents/imi_projects/MBV/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.MR_DWI_reg.nii.gz'),
            'DWI')

        self.content.add_tab(to_input=False, label='test', image=sitk.ReadImage(
            '/home/paul/Documents/imi_projects/MBV/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.MR_T1_reg.nii.gz'))

        self.add_image(sitk.ReadImage(
            '/home/paul/Documents/imi_projects/MBV/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.MR_T1_reg.nii.gz'),
            'T1')

    def init_menuBar(self):
        menubar = self.menuBar()

        exit_action = QAction('Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close)

        open_action = QAction('Load image', self)
        open_action.setShortcut('Ctrl+L')
        open_action.setStatusTip('Load MR image into currently selected input widget')
        open_action.triggered.connect(self.handler.load_image)

        open_mask_action = QAction('Load reference mask', self)
        open_mask_action.setShortcut('Ctrl+M')
        open_mask_action.setStatusTip('Load reference mask of the image')

        file_menu = menubar.addMenu('&File')
        file_menu.addActions([open_action, open_mask_action, exit_action])

        view_menu = menubar.addMenu('&View')
        orientation_submenu = view_menu.addMenu('Change image &orientation')
        xy_orientation_action = QAction('XY plane', self)
        xy_orientation_action.setStatusTip('Set orientation to XY (axial) plane')
        xz_orientation_action = QAction('XZ plane', self)
        xz_orientation_action.setStatusTip('Set orientation to XZ (coronal) plane')
        yz_orientation_action = QAction('YZ plane', self)
        yz_orientation_action.setStatusTip('Set orientation to YZ (sagittal) plane')

        orientation_submenu.addActions([xy_orientation_action, xz_orientation_action, yz_orientation_action])

        marker_submenu = view_menu.addMenu('&Marker')
        clear_markers_action = QAction('&Clear all markers', self)
        clear_markers_action.setStatusTip('Clear all markers from all images')
        remove_last_marker_action = QAction('Remove &last marker', self)
        remove_last_marker_action.setStatusTip('Remove the last placed marker')

        marker_submenu.addActions([clear_markers_action, remove_last_marker_action])

    def add_image(self, image, name):
        self.content.add_input_tab(name, image)


class ProjectMainWindowContent(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(ProjectMainWindowContent, self).__init__(parent)

        self.input_tabs = ImageTabWidget(self)
        self.output_tabs = ImageTabWidget(self)

        central_hbox = QtWidgets.QHBoxLayout(self)
        central_hbox.addWidget(self.input_tabs)
        central_hbox.addWidget(self.output_tabs)
        central_hbox.addWidget(ControlPanel(self))
        self.setLayout(central_hbox)

    def add_input_tab(self, label, image, at_index=None):
        self.add_tab(to_input=True, label=label, image=image, at_index=at_index)

    def add_output_tab(self, label, image, at_index=None):
        self.add_tab(to_input=False, label=label, image=image, at_index=at_index)

    def add_tab(self, to_input, label, image, at_index=None):
        tab_widget = self.input_tabs if to_input else self.output_tabs

        viewer = ImageViewer(parent=self.input_tabs, itk_image=image)
        if at_index is None:
            tab_index = tab_widget.addTab(viewer, label)
        else:
            tab_index = tab_widget.insertTab(at_index, viewer, label)

        return tab_index


class ImageTabWidget(QtWidgets.QTabWidget):
    def __init__(self, parent=None):
        super(ImageTabWidget, self).__init__(parent)
        self.setTabsClosable(True)
        self.tabCloseRequested.connect(self.removeTab)
        self.currentChanged.connect(self.hide_other_widgets)

    def removeTab(self, index: int) -> None:
        # override the close button method with the deletion of the contained widget
        self.widget(index).deleteLater()

    def hide_other_widgets(self, index):
        cur = 0
        while self.widget(cur):
            if cur != index:
                self.widget(cur).setVisible(False)
            else:
                self.widget(cur).setVisible(True)

            cur += 1




# class ImageWidget(QtWidgets.QWidget):
#     def __init__(self, parent=None):
#         super(ImageWidget, self).__init__(parent)
#
#         self.frame = QtWidgets.QFrame()
#
#         self.vl = QtWidgets.QVBoxLayout()
#         self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
#         self.vl.addWidget(self.vtkWidget)
#
#         self.slice_slider = QtWidgets.QSlider(orientation=1, parent=self)
#         self.vl.addWidget(self.slice_slider)
#
#         self.ren = vtk.vtkRenderer()
#         self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
#         self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
#
#         # Create source
#         source = vtk.vtkSphereSource()
#         source.SetCenter(0, 0, 0)
#         source.SetRadius(5.0)
#
#         # Create a mapper
#         mapper = vtk.vtkPolyDataMapper()
#         mapper.SetInputConnection(source.GetOutputPort())
#
#         # Create an actor
#         actor = vtk.vtkActor()
#         actor.SetMapper(mapper)
#
#         self.ren.AddActor(actor)
#
#         self.ren.ResetCamera()
#
#         self.frame.setLayout(self.vl)
#
#         self.iren.Initialize()


class ControlPanel(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(ControlPanel, self).__init__(parent)

        self.content_vbox = QtWidgets.QVBoxLayout()

        # add 6 LineEdits for value inputs
        self.line_edits = [QtWidgets.QLineEdit() for _ in range(6)]
        self.init_valueinputs()

        # add buttons
        self.button_preprocessing = QtWidgets.QPushButton('Pre-Processing', self)
        self.button_segmentation = QtWidgets.QPushButton('Segmentation', self)
        self.button_postprocessing = QtWidgets.QPushButton('Post-Processing', self)
        self.button_evaluation = QtWidgets.QPushButton('Evaluation', self)

        self.button1 = get_button_bound_to_function('Button 1', lambda *args: None, self)
        # self.button1.clicked.connect(self.print_le1)
        self.button2 = QtWidgets.QPushButton('Button 2', self)
        self.button3 = QtWidgets.QPushButton('Button 3', self)
        self.button4 = QtWidgets.QPushButton('Button 4', self)
        self.button5 = QtWidgets.QPushButton('Button 5', self)
        self.button6 = QtWidgets.QPushButton('Button 6', self)

        # TODO: figure out the binding of buttons
        self.init_buttons()

        self.content_vbox.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.setLayout(self.content_vbox)

    def init_valueinputs(self):
        line_edit_form = QtWidgets.QFormLayout(self)

        for i, le in enumerate(self.line_edits):
            le.setValidator(QDoubleValidator())
            #le.textChanged.connect(self.print_lineedit_value)
            line_edit_form.addRow('Value {}'.format(i+1), le)

        widget = QtWidgets.QWidget(self)
        widget.setLayout(line_edit_form)
        self.content_vbox.addWidget(widget)

    def init_buttons(self):
        for c in self.findChildren(QtWidgets.QPushButton):
            self.content_vbox.addWidget(c)

    def get_value(self, index: int):
        assert 0 < index <= len(self.line_edits), 'Tried to get value number {}. You need to choose one between 1 and {}.'.format(index, len(self.line_edits))
        # replace comma with point in text before converting to float
        text_value = self.line_edits[index - 1].text().replace(',', '.')

        try:
            float_value = float(text_value)
        except ValueError:
            # on an empty string, return 0
            float_value = 0

        return float_value


def get_button_bound_to_function(btn_label, func, btn_parent=None):
    btn = QtWidgets.QPushButton(btn_label, btn_parent)
    btn.clicked.connect(func)
    return btn


if __name__ == '__main__':
    pass
