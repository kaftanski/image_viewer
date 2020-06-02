import SimpleITK as sitk
from math import floor
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


COLORS = {
    'red': [1, 0, 0],
    'green': [0, 1, 0],
    'blue': [0, 0, 1]
}


class ImageViewer(QtWidgets.QWidget):
    # TODO:
    #  image
    #  mask
    #  markers
    #  slice-slider (if possible. else in QWidget)
    #  different orientations of slices
    def __init__(self, itk_image=None, parent=None):
        super(ImageViewer, self).__init__(parent)
        self.viewportsize_x = 800
        self.viewportsize_y = 800
        self.masks = []  # maybe dict todo
        self.markers = []

        # setup the image viewer

        # first the Qt part as an embedding
        # self.image_frame = QtWidgets.QFrame(parent=self if parent is None else parent)
        self.vtk_widget = QVTKRenderWindowInteractor(parent=self)

        # now the vtk part for the actual image rendering
        # test_reader = vtk.vtkNIFTIImageReader()
        # test_reader.SetFileName('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.MR_DWI_reg.nii.gz')
        # test_reader.Update()
        # self.image = test_reader.GetOutput()

        if itk_image is None:
            # load example image
            itk_image = sitk.ReadImage('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.MR_Flair_reg.nii.gz')

        self.itk_image = itk_image
        self.itk_to_vtk_image_filter = ITKtoVTKImageFilter()
        self.image = self.itk_to_vtk_image_filter(itk_image)

        # greyval_range = self.image.GetScalarRange()
        # self.wl_lut = vtk.vtkWindowLevelLookupTable()
        # self.wl_lut.SetWindow(greyval_range[1] - greyval_range[0])
        # self.wl_lut.SetLevel((greyval_range[1] - greyval_range[0]) / 2)
        # self.wl_lut.Build()

        # self.ds_mapper = vtk.vtkDataSetMapper()
        # self.ds_mapper.SetInputDataObject(self.image)
        # self.ds_mapper.SetLookupTable(self.wl_lut)
        # self.ds_mapper.SetScalarRange(0, 255)

        # self.image_actor = vtk.vtkActor()
        # self.image_actor.SetMapper(self.ds_mapper)

        self.blender = vtk.vtkImageBlend()
        # self.add_mask_image(None)
        self.blender.AddInputDataObject(self.image)
        # self.blender.AddInputConnection(self.ds_mapper.GetOutputPort())
        self.blender.SetBlendModeToCompound()

        renwin = self.vtk_widget.GetRenderWindow()

        self.vtk_image_viewer = vtk.vtkImageViewer2()
        self.vtk_image_viewer.SetSize(self.viewportsize_x, self.viewportsize_y)
        self.vtk_image_viewer.SetupInteractor(self.vtk_widget)
        self.vtk_image_viewer.SetRenderWindow(renwin)

        self.renwin_interactor = self.vtk_widget.GetRenderWindow().GetInteractor()
        interactor_style = IMIProjectInteractorStyleImage(self)
        self.renwin_interactor.SetInteractorStyle(interactor_style)
        interactor_style.SetInteractor(self.renwin_interactor)

        # self.vtk_image_viewer.SetInputConnection(self.blender.GetOutputPort())
        self.vtk_image_viewer.SetInputData(self.image)
        self.vtk_image_viewer.SetSlice(60)
        greyval_range = self.image.GetScalarRange()
        self.vtk_image_viewer.SetColorWindow(greyval_range[1] - greyval_range[0])
        self.vtk_image_viewer.SetColorLevel((greyval_range[1] - greyval_range[0]) / 2)
        self.vtk_image_viewer.Render()

        # define a slider for image slicing
        self.slice_slider = QtWidgets.QSlider(orientation=Qt.Vertical, parent=self)
        self.init_slider()

        # define a text field to display image coordinates
        coordinates_hbox = QtWidgets.QHBoxLayout()
        self.pixel_info_label = PixelInfoQLabel(parent=self)
        coordinates_hbox.addWidget(self.pixel_info_label)

        # layout the QtWidget
        main_layout = QtWidgets.QVBoxLayout()
        image_with_slider_layout = QtWidgets.QHBoxLayout()
        image_with_slider_layout.addWidget(self.vtk_widget)
        image_with_slider_layout.addWidget(self.slice_slider)
        main_layout.addLayout(image_with_slider_layout)
        main_layout.addLayout(coordinates_hbox)
        self.setLayout(main_layout)

        # self.add_mask_image(None)

        self.vtk_widget.Start()
        self.show()
        self.renwin_interactor.Initialize()

    def update_slice(self, new_slice):
        self.vtk_image_viewer.SetSlice(new_slice)
        self.pixel_info_label.set_coordinate(self.vtk_image_viewer.GetSliceOrientation(), new_slice)
        self.show_pixel_info(self.pixel_info_label.coords)
        self.update_masks()

    def init_slider(self):
        self.slice_slider.setRange(self.vtk_image_viewer.GetSliceMin(), self.vtk_image_viewer.GetSliceMax())
        self.slice_slider.setValue(self.vtk_image_viewer.GetSlice())
        self.slice_slider.valueChanged.connect(self.update_slice)

    def move_slice_forward(self):  # TODO: change this ugly workaround
        # moves the slider, which in turn changes the slice
        self.slice_slider.setValue(self.vtk_image_viewer.GetSlice() + 1)
        # self.update_slice(self.vtk_image_viewer.GetSlice() + 1)

    def move_slice_backward(self):
        # moves the slider, which in turn changes the slice
        self.slice_slider.setValue(self.vtk_image_viewer.GetSlice() - 1)
        # self.update_slice(self.vtk_image_viewer.GetSlice() - 1)

    def add_marker(self, position):
        assert len(position) == 3
        marker = ImageMarker(position)
        self.markers.append(marker)

        ren = self.vtk_widget.GetRenderWindow().GetRenderers().GetFirstRenderer()
        ren.AddActor(marker.GetActor())
        ren.ResetCamera()
        # TODO remember markerlist

    def add_mask_image(self, vtk_binary_image):
        # TODO assert spacing and dimensions etc equal
        #assert vtk_binary_image.GetExtent() == self.vtk_image_viewer.GetInput().GetExtent()

        test_mask = ImageMask(None, 0.3, 'blue', mask_scalar_type=self.image.GetScalarType())
        self.masks.append(test_mask)

        label_overlay = sitk.LabelOverlayImageFilter()
        overlayed = label_overlay.Execute(self.itk_image, sitk.ReadImage('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.OT_reg.nii.gz'))
        overlayed = self.itk_to_vtk_image_filter(overlayed)
        self.vtk_image_viewer.SetInputDataObject(overlayed)

        # self.vtk_image_viewer.GetRenderer().AddActor(test_mask.GetActor())
        # self.blender.AddInputConnection(test_mask.GetOutputPort())

        # self.update_masks()

    def update_masks(self):
        for i, mask in enumerate(self.masks):
            mask_actor = mask.GetActor()
            mask_actor.SetDisplayExtent(self.vtk_image_viewer.GetImageActor().GetDisplayExtent())

            self.blender.ReplaceNthInputConnection(i+1, mask.GetOutputPort())
            self.blender.Update()
            # self.vtk_image_viewer.Render()
            #if not self.vtk_image_viewer.GetRenderer().GetActors().IsItemPresent(mask_actor):
            # print(self.vtk_image_viewer.GetImageActor())
            # print(self.vtk_image_viewer.GetRenderer().GetActors())

            # self.vtk_image_viewer.GetRenderer().AddActor(mask_actor)
            # self.vtk_image_viewer.Render()
            # interactor_style.GetDefaultRenderer().ResetCamera()

    def set_image(self, image):
        self.image = image
        self.init_slider()

    def show_pixel_info(self, coords):
        if coords is None:
            self.pixel_info_label.clear()
        else:
            x, y, z = [floor(ind) for ind in coords]
            self.pixel_info_label.set_values(x, y, z, self.image.GetScalarComponentAsFloat(x, y, z, 0))


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


class ImageMask:
    def __init__(self, itk_binary_image, alpha=0.3, color='red', mask_scalar_type=5):  # default type is 5 = ushort
        self.itk_mask = itk_binary_image
        self.alpha = alpha
        self.color = COLORS[color]

        self.lookup_table = vtk.vtkLookupTable()
        self.lookup_table.SetNumberOfTableValues(2)
        self.lookup_table.SetRange(0.0, 1.0)
        self.lookup_table.SetTableValue(0, 0.0, 0.0, 0.0, 0.0)  # label 0 is transparent
        rgba = [*self.color, alpha]
        self.lookup_table.SetTableValue(1, rgba)  # label 1 is specified by color and alpha values
        self.lookup_table.Build()

        self.label_map = vtk.vtkImageMapToColors()
        self.label_map.SetLookupTable(self.lookup_table)
        self.label_map.PassAlphaToOutputOn()

        test_reader = vtk.vtkNIFTIImageReader()
        test_reader.SetFileName('/home/paul/Documents/imi_projects/MBV/Projekt/MIPImages/ISLES2015_Train/01/VSD.Brain.01.O.OT_reg.nii.gz')
        test_reader.Update()

        self.label_map.SetInputData(test_reader.GetOutput())
        # label_map.Update()

        self.image_cast = vtk.vtkImageCast()
        self.image_cast.SetOutputScalarType(mask_scalar_type)
        self.image_cast.SetInputConnection(self.label_map.GetOutputPort())
        self.image_cast.Update()

        self.mask_image_actor = vtk.vtkImageActor()
        self.mask_image_actor.GetMapper().SetInputData(self.image_cast.GetOutput())
        # print(self.mask_image_actor.GetMapper().GetInput())
        self.mask_image_actor.Update()

    def GetActor(self):
        return self.mask_image_actor

    def GetOutputPort(self):
        return self.mask_image_actor.GetMapper().GetOutputPort()
        # return self.image_cast.GetOutputPort()

    def GetOutputDataObject(self):
        result = self.mask_image_actor.GetMapper().GetOutputDataObject(0)
        return result


class ITKtoVTKImageFilter:
    """ convert sitk image to vtk image via the raw data in a numpy array """
    def __init__(self):
        self.input_image = None
        self.output_image = None
        self.vtk_image_importer = vtk.vtkImageImport()
        self.image_as_array = None

    def __call__(self, input_image):
        self.input_image = input_image
        self.Execute()
        return self.output_image

    def Execute(self):
        self.image_as_array = sitk.GetArrayFromImage(self.input_image)
        # resulting array is indexed as (z, y, x), which has to be reversed before importing
        image_as_string = self.image_as_array.flatten().tostring()
        self.vtk_image_importer.CopyImportVoidPointer(image_as_string, len(image_as_string))
        self.vtk_image_importer.SetDataScalarTypeToUnsignedShort()  # todo
        self.vtk_image_importer.SetNumberOfScalarComponents(1)  # image only contains intensity values in one channel
        # todo: generalize to multiple channels (rgb/rgba)?
        # self.vtk_image_importer.SetDataDirection(self.input_image.GetDirection())  # TODO: no data direction in wrapped vtk?

        # set the image extent with reversed axes compared to the numpy array
        self.vtk_image_importer.SetWholeExtent([self.image_as_array.shape[-i // 2] - 1 if i % 2 else 0
                                                for i in range(2 * len(self.image_as_array.shape))])
        self.vtk_image_importer.SetDataExtentToWholeExtent()

        self.vtk_image_importer.SetDataSpacing(self.input_image.GetSpacing())
        self.vtk_image_importer.SetDataOrigin(self.input_image.GetOrigin())

        self.vtk_image_importer.Update()

        self.output_image = self.vtk_image_importer.GetOutput()

    def GetOutput(self):
        if self.output_image is None:
            raise UserWarning('Tried to get output of ITKVTKImageFilter, which has not been executed.')

        return self.output_image


# TODO: check dimensions!!
# class CanvasImageMarker:
#     def __init__(self, world_position, extent, radius=10):
#         assert len(world_position) == 3
#         self.world_position = world_position
#
#         self.drawing = vtk.vtkImageCanvasSource2D()
#         self.drawing.SetExtent(extent)
#         self.drawing.SetDrawColor(255, 255, 255)
#
#         self.drawing.DrawCircle(round(world_position[0]), round(world_position[1]), radius)
#         # self.drawing.DrawPoint(round(world_position[0]), round(world_position[1]))


class ImageMarker:
    def __init__(self, world_position, radius=3):
        assert len(world_position) == 3
        self.world_position = world_position

        sphere_source = vtk.vtkSphereSource()
        sphere_source.SetRadius(radius)
        sphere_source.SetThetaResolution(50)
        sphere_source.SetPhiResolution(50)
        sphere_source.SetCenter(self.world_position)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere_source.GetOutputPort())

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().SetColor(COLORS['red'])
        self.actor.GetProperty().SetAmbient(0.7)
        self.actor.GetProperty().SetDiffuse(0.0)
        self.actor.GetProperty().SetSpecularPower(5.0)
        self.actor.GetProperty().SetSpecular(0.0)

    def GetActor(self):
        return self.actor


class IMIProjectInteractorStyleImage(vtk.vtkInteractorStyleImage):
    def __init__(self, imi_image_viewer):
        self.imi_image_viewer = imi_image_viewer
        self.vtk_image_viewer = imi_image_viewer.vtk_image_viewer
        self.SetDefaultRenderer(self.vtk_image_viewer.GetRenderer())

        self.AddObserver('MouseWheelForwardEvent', self.on_mouse_wheel_forward)  # TODO: remove observers?
        self.AddObserver('MouseWheelBackwardEvent', self.on_mouse_wheel_backward)
        self.AddObserver('RightButtonPressEvent', self.AddMarker)
        self.AddObserver('MouseMoveEvent', self.on_mouse_move)

        self.key_to_function_mapping = {
            'Up': self.on_mouse_wheel_forward,
            'Down': self.on_mouse_wheel_backward,
            'x': (lambda *args, **kwargs: None)
        }
        self.AddObserver('KeyPressEvent', self.HandleKeyPress, -1.0)

    def on_mouse_wheel_forward(self, *args, **kwargs):
        self.imi_image_viewer.move_slice_forward()
        # self.imi_image_viewer.show_pixel_info(self.get_event_image_coordinates())

    def on_mouse_wheel_backward(self, *args, **kwargs):
        self.imi_image_viewer.move_slice_backward()
        # self.imi_image_viewer.show_pixel_info(self.get_event_image_coordinates())

    def on_mouse_move(self, *args, **kwargs):
        position = self.get_event_image_coordinates()
        self.imi_image_viewer.show_pixel_info(position)

        self.OnMouseMove()  # cannot be overridden, but only used by subclass in python

    def get_event_image_coordinates(self):
        event_x, event_y = self.GetInteractor().GetEventPosition()

        # print('Transforming event at ({}, {}) into image world coords'.format(event_x, event_y))
        dimensions = self.imi_image_viewer.image.GetDimensions()

        # pick the
        picker = vtk.vtkCellPicker()
        picker.Pick(event_x, event_y, 0, self.GetDefaultRenderer())
        pixel_position = list(picker.GetPickPosition())

        # don't do anything if clicked out of image bounds
        if any(pixel_position[i] < 0 or pixel_position[i] >= dimensions[i] for i in range(2)):
            # print('Event out of image bounds')
            return None

        # orientation is defined as one of the following:
        #     vtk.vtkImageViewer2.SLICE_ORIENTATION_XY = 2,
        #     vtk.vtkImageViewer2.SLICE_ORIENTATION_XZ = 1,
        #     vtk.vtkImageViewer2.SLICE_ORIENTATION_YZ = 0
        # depending on the current view.
        orientation = self.vtk_image_viewer.GetSliceOrientation()

        # assign the slice index to the right dimension given as orientation
        pixel_position[orientation] = self.vtk_image_viewer.GetSlice() * self.imi_image_viewer.image.GetSpacing()[
            orientation] + self.imi_image_viewer.image.GetOrigin()[orientation]  # TODO: warum origin?

        return pixel_position

    def HandleKeyPress(self, *args, **kwargs):
        key = self.GetInteractor().GetKeySym()
        print('pressed key {}'.format(key))

        try:
            self.key_to_function_mapping[key](args, kwargs)
        except KeyError:
            self.OnKeyPress()

    def AddMarker(self, *args, **kwargs):
        image_coords = self.get_event_image_coordinates()

        if image_coords is not None:
            print('Placing marker at coordinates:', image_coords)
            marker = ImageMarker(image_coords)

            renderer = self.GetDefaultRenderer()
            renderer.AddActor(marker.GetActor())
            renderer.ResetCamera()
            self.vtk_image_viewer.Render()  # TODO: keep zoom of image

            # canvas = CanvasImageMarker(image_coords, self.image_viewer.GetInput().GetExtent())
            #
            # canvas.drawing.SetScalarType(self.image_viewer.GetInput().GetScalarType())
            # blender = vtk.vtkImageBlend()
            # image_as_canvas = vtk.vtkImageCanvasSource2D()
            # image_as_canvas.InitializeCanvasVolume(self.image_viewer.GetInput())
            # blender.AddInputDataObject(0, self.image_viewer.GetInput())  # port 0
            # blender.AddInputConnection(1, canvas.drawing.GetOutputPort())  # port 1
            # blender.SetOpacity(0, 0.7)
            # blender.SetOpacity(1, 0.3)
            # self.image_viewer.SetInputConnection(blender.GetOutputPort())
            # self.image_viewer.Render()

        self.OnRightButtonDown()


if __name__ == '__main__':
    import sys

    app = QtWidgets.QApplication(sys.argv)
    w = ImageViewer()
    sys.exit(app.exec_())
