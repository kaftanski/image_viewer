import SimpleITK as sitk
import numpy as np
from matplotlib.colors import ListedColormap


def get_3d_plane_index(slice_index, orientation):
    index = [slice(0, None), slice(0, None)]  # get the whole dimensions of planar axes
    index.insert(orientation, slice_index)
    return tuple(index)


def get_aspect_ratio_for_plane(spacing, orientation, image_dimensions):
    dims = list(d for d in range(3) if d != orientation)
    ratio = spacing[dims[1]] * image_dimensions[1] / (spacing[dims[0]] * image_dimensions[0])
    return ratio


def compatible_metadata(image1: sitk.Image, image2: sitk.Image, check_size=True, check_spacing=True, check_origin=True):
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


def index_compatibility(index):
    # change to the other method (sitk uses x, y, z, numpy uses z, y, x
    assert len(index) == 3
    return index[-1::-1]


def pix2world(index, spacing):
    assert len(index) == len(spacing)
    return tuple(i*s for i, s in zip(index, spacing))


def world2pix(index, spacing):
    assert len(index) == len(spacing)
    return tuple(i/s for i, s in zip(index, spacing))


def add_mask_to_image(ax, mask, aspect, alpha=0.3, color='r'):
    """ can be multilabel, if mask is given as (x, y, l) with l different binary images """
    if mask.sum() == 0:
        return None

    cmap = ListedColormap([color])
    mask = np.ma.masked_where(mask == 0, mask)
    plot = ax.matshow(mask, aspect=aspect, cmap=cmap, alpha=alpha)

    return plot
