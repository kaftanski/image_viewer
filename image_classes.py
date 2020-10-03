from typing import List, Union, Sequence, Tuple

import SimpleITK as sitk
import numpy as np

from image_utils import index_compatibility, get_3d_plane_index


class Image:
    def __init__(self, image_array: np.ndarray, spacing: Sequence[float] = None, origin: Sequence[Union[int, float]] = None):
        # ensure 2D or 3D image
        image_dim = len(image_array.shape)
        assert 2 <= image_dim <= 3

        self.image_array = image_array
        self.spacing = spacing if spacing is not None else np.ones(image_dim)
        self.origin = origin if origin is not None else np.zeros(image_dim)

    def GetImageArray(self):
        """
        @return: the numpy array containing the image intensities
        """
        return self.image_array

    def GetSpacing(self) -> Sequence[float]:
        """
        @return: the image spacing
        """
        return self.spacing

    def GetOrigin(self) -> Sequence[Union[int, float]]:
        """
        @return: the image origin
        """
        return self.origin

    def GetSize(self) -> Tuple[int]:
        """
        @return: the dimensions of the image
        """
        return self.image_array.shape


class ImageMask:
    """ image mask container holding a mask image and display parameters """
    def __init__(self, sitk_binary_image: sitk.Image, alpha: float = 0.3, color: str = 'r'):
        """

        @param sitk_binary_image: the image interpreted as a mask
        @param alpha: the alpha value for the mask overlay (between 0 and 1)
        @param color: the desired color for the mask (all matplotlib's color definitions are possible)
        """
        self.itk_mask = sitk_binary_image
        self.mask_as_array = sitk.GetArrayFromImage(self.itk_mask)
        self.alpha = alpha
        self.color = color

    def get_slice(self, slice_index: int, orientation: int) -> np.ndarray:
        """ Get a slice from the mask

        @param slice_index: the index in the slicing dimension
        @param orientation: the slicing dimension
        @return: a slice from the dimension given in orientation at the given index
        """
        return self.mask_as_array[index_compatibility(get_3d_plane_index(slice_index, orientation))]

    def get_spacing(self) -> List[float]:
        """

        @return: the spacing of the mask image
        """
        return self.itk_mask.GetSpacing()


class ImageMarker:
    """ image marker container holding a pixel position and the default parameters for displaying markers """
    # TODO: kind of an unnecessary class

    STANDARD_COLOR = 'b'  # follows matplotlib colors definitions
    STANDARD_SIZE = 10  # size of markers in pt

    def __init__(self, pixel_position: Sequence[Union[int, float]]):
        """
        @param pixel_position: 3-dimensional (possibly continouus) index of the marker
        """
        assert len(pixel_position) == 3
        self.pixel_position = pixel_position

    def __str__(self):
        """
        @return: a nicely formatted string displaying the marker's coordinates
        """
        return 'ImageMarker at pixel position ({}, {}, {})'.format(*self.pixel_position)
