import SimpleITK as sitk

from image_utils import index_compatibility, get_3d_plane_index


class ImageMask:
    def __init__(self, itk_binary_image, alpha=0.3, color='r'):
        self.itk_mask = itk_binary_image
        self.mask_as_array = sitk.GetArrayFromImage(self.itk_mask)
        self.alpha = alpha
        self.color = color

    def get_slice(self, slice_index, orientation):
        return self.mask_as_array[index_compatibility(get_3d_plane_index(slice_index, orientation))]

    def get_spacing(self):
        return self.itk_mask.GetSpacing()


class ImageMarker:
    # TODO: kind of an unnecessary class
    STANDARD_COLOR = 'b'  # follows matplotlib colors definitions
    STANDARD_SIZE = 10  # size of markers in pt

    def __init__(self, pixel_position):
        assert len(pixel_position) == 3
        self.pixel_position = pixel_position

    def __str__(self):
        return 'ImageMarker at pixel position ({}, {}, {})'.format(*self.pixel_position)
