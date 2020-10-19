from typing import Tuple, Union, Sequence

import SimpleITK as sitk
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap
from matplotlib.image import AxesImage


class Image:
    def __init__(self, image_array: np.ndarray, spacing: Sequence[float] = None, origin: Sequence[Union[int, float]] = None, reverse_indexing: bool = False):
        # ensure 2D or 3D image
        image_dim = len(image_array.shape)
        assert 2 <= image_dim <= 3

        self.image_array = image_array
        self.spacing = spacing if spacing is not None else np.ones(image_dim)
        self.origin = origin if origin is not None else np.zeros(image_dim)
        self.reverse_indexing = reverse_indexing

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
        size = list(self.image_array.shape)
        if self.reverse_indexing:
            size = size[::-1]
        return tuple(size)

    def GetPixel(self, x: int, y: int, z: int = None) -> Union[int, float]:
        if z is None:
            assert len(self.GetSize()) == 2, 'Error: Trying to access a 3D image with only 2 indices. ' \
                                             'If needed, use numpy-like indexing with []-brackets.'
            return self.image_array[y, x] if self.reverse_indexing else self.image_array[x, y]
        else:
            assert len(self.GetSize()) == 3, 'Error: Trying to access a 2D image with 3 indices.'
            return self.image_array[z, y, x] if self.reverse_indexing else self.image_array[x, y, z]

    def SetPixel(self, value: Union[int, float], x: int, y: int, z: int = None):
        if z is None:
            assert len(self.GetSize()) == 2, 'Error: Trying to access a 3D image with only 2 indices. ' \
                                             'If needed, use numpy-like indexing with []-brackets.'

            if self.reverse_indexing:
                self.image_array[y, x] = value
            else:
                self.image_array[x, y] = value
        else:
            assert len(self.GetSize()) == 3, 'Error: Trying to access a 2D image with 3 indices.'

            if self.reverse_indexing:
                self.image_array[z, y, x] = value
            else:
                self.image_array[x, y, z] = value

    def __getitem__(self, item) -> np.ndarray:
        index_compatibility(item)
        return self.image_array[item]

    def __setitem__(self, key, value):
        self.image_array[key] = value


class ImageMask:
    """ image mask container holding a mask image and display parameters """
    def __init__(self, binary_image: sitk.Image, alpha: float = 0.3, color: str = 'r'):
        """

        @param binary_image: the image interpreted as a mask
        @param alpha: the alpha value for the mask overlay (between 0 and 1)
        @param color: the desired color for the mask (all matplotlib's color definitions are possible)
        """
        self.mask_image = binary_image
        self.mask_as_array = sitk.GetArrayFromImage(self.mask_image)
        self.alpha = alpha
        self.color = color

    def get_slice(self, slice_index: int, orientation: int) -> np.ndarray:
        """ Get a slice from the mask

        @param slice_index: the index in the slicing dimension
        @param orientation: the slicing dimension
        @return: a slice from the dimension given in orientation at the given index
        """
        return self.mask_as_array[index_compatibility(get_3d_plane_index(slice_index, orientation))]

    def get_spacing(self) -> Sequence[float]:
        """

        @return: the spacing of the mask image
        """
        return self.mask_image.GetSpacing()


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


def get_3d_plane_index(slice_index: int, orientation: int) -> Tuple[Union[slice, int]]:
    """ Returns an index to extract a certain slice from the given dimension in orientation

    @param slice_index: the index of the desired slice
    @param orientation: the sliceing dimension
    @return: an index to get one slice
    """
    index = [slice(0, None), slice(0, None)]  # get the whole dimensions of planar axes
    index.insert(orientation, slice_index)
    return tuple(index)


def get_aspect_ratio_for_plane(spacing: Sequence[float], orientation: int, image_dimensions: Sequence[int]) -> float:
    """ Computes the necessary aspect ratio for two axes given a spacing

    @param spacing: the spacing for all 3 dimensions
    @param orientation: the current slicing dimension
    @param image_dimensions: the size of all dimensions
    @return: the aspect ratio between the two plane dimensions
    """
    dims = list(d for d in range(3) if d != orientation)
    ratio = spacing[dims[1]] * image_dimensions[1] / (spacing[dims[0]] * image_dimensions[0])
    return ratio


def compatible_metadata(image1: Image, image2: Image, check_size: bool = True, check_spacing: bool = True, check_origin: bool = True) -> bool:
    """ Compares the metadata of two images and determines if all checks are successful or not.
    Comparisons are carried out with a small tolerance (0.0001).

    @param image1: first image
    @param image2: second image
    @param check_size: if true, check if the sizes of the images are equal
    @param check_spacing: if true, check if the spacing of the images are equal
    @param check_origin: if true, check if the origin of the images are equal
    @return: true, if images are equal in all given checks, false if one of them failed
    """
    all_parameters_equal = True
    tolerance = 1e-4

    if check_size:
        size1 = image1.GetSize()
        size2 = image2.GetSize()
        if size1 != size2:
            all_parameters_equal = False
            print(f'Images do not have the same size ({size1} != {size2})')

    if check_spacing:
        spacing1 = image1.GetSpacing()
        spacing2 = image2.GetSpacing()
        if any(list(abs(s1-s2) > tolerance for s1, s2 in zip(spacing1, spacing2))):
            all_parameters_equal = False
            print(f'Images do not have the same spacing ({spacing1} != {spacing2})')

    if check_origin:
        origin1 = image1.GetOrigin()
        origin2 = image2.GetOrigin()
        if any(list(abs(o1-o2) > tolerance for o1, o2 in zip(origin1, origin2))):
            all_parameters_equal = False
            print(f'Images do not have the same origin ({origin1} != {origin2})')

    return all_parameters_equal


def index_compatibility(index: Sequence[int]) -> Sequence[int]:
    """ Convert an index for numpy into one for SimpleITK and vice-versa.
    Background: numpy uses (z, y, x) indexing, SimpleITK uses (x,y,z).
    This function therefore reverses the index.

    @param index: the index to be adapted for the other indexing scheme
    @return: the index for the other indexing scheme
    """
    # change to the other method (sitk uses x, y, z, numpy uses z, y, x)
    assert len(index) == 3
    return index[::-1]


def pix2world(pixel_coords: Sequence[Union[int, float]], spacing: Sequence[float]) -> Tuple[float]:
    """ Converts pixel coordinates into world coordinates using the given voxel spacing.

    @param pixel_coords: the pixel coordinates
    @param spacing: the voxel spacing of an image
    @return: the continuous world coordinates
    """
    assert len(pixel_coords) == len(spacing)
    return tuple(i * s for i, s in zip(pixel_coords, spacing))


def world2pix(world_coords: Sequence[Union[int, float]], spacing: Sequence[float]) -> Tuple[float]:
    """ Converts world coordinates into pixel coordinates using the given voxel spacing.

    @param world_coords: the world coordinates
    @param spacing: the voxel spacing of an image
    @return: the continuous pixel coordinates
    """
    assert len(world_coords) == len(spacing)
    return tuple(i / s for i, s in zip(world_coords, spacing))


def add_mask_to_image(ax: Axes, mask: np.ndarray, aspect: float, alpha: float = 0.3, color: str = 'r') -> Union[AxesImage, None]:
    """ Add a single label mask an an alpha channel to the given axis.
    Any value in the mask not equal to 0 is considered as object, pixels with value 0 are considered as background.

    @param ax: the matplotlib Axes where the mask will be shown
    @param mask: a 2d array containing the mask image.
                 0-pixels are considered background, other pixels are drawn as foreground
    @param aspect: the aspect ratio of the axes
    @param alpha: the alpha value for the mask
    @param color: the color of the mask (use a string matplotlib.colors can interpret)
    @return: the AxesImage resulting from the plot call
    """
    if mask.sum() == 0:
        return None

    cmap = ListedColormap([color])
    mask = np.ma.masked_where(mask == 0, mask)
    plot = ax.matshow(mask, aspect=aspect, cmap=cmap, alpha=alpha)

    return plot
