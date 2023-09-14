# 3D Slice Image Viewer
Interactive image viewer for 3D medical images that can be used with inline functions (``image_viewing.py``) or through a command line interface (``viewImage.py``).
The viewer is a PyQt5 application with a matplotlib canvas to draw the images.
Pyplot-based image viewer for 3D medical images that was developed for a Bachelor's project in Medical Image Processing at the Institute of Medical Informatics (University of Lübeck).

## Image viewing
Open an image viewer with ``image_viewing.show_image``. You can scroll through slices across the three principal axes of the 3D image with the mouse.

![image_viewer](https://user-images.githubusercontent.com/49405884/220111359-ceac8ea9-3917-49c7-945b-9ee7d8fafcda.png)

## View segmentation masks
Open an image with corresponding mask image (different integer values are colored differently): ``image_viewing.show_image_with_mask``

![image_segmentation_example](https://user-images.githubusercontent.com/49405884/220112222-b313dcab-05cb-431a-bc48-0f226eec63c4.png)

## Interactive Seed Points
Set markers or seed points interactively and return them on closing the window (``image_viewing.show_and_return_markers``). This is useful for interactive region growing pipelines.

![image_viewer_markers](https://user-images.githubusercontent.com/49405884/220112488-6baed0f3-c4d0-4ffc-be7a-af9be9afe685.png)

## Requirements
* Python3
* matplotlib
* SimpleITK
* PyQt5
