import SimpleITK as sitk
import sys
import image_viewing as vis

assert len(sys.argv) > 1, 'No input image specified'

# load example image from argv
input_image = sitk.ReadImage(sys.argv[1])
vis.show_image(input_image, 'Input Image (non-blocking window)', blocking=False)

# pre-processing
print('smoothing')
smoothed = sitk.DiscreteGaussian(input_image, 2)
vis.show_image(smoothed, 'Smoothed Image', blocking=False)

# get seedpoints
seeds = vis.show_and_return_markers(smoothed, 'Set Seedpoints')
print(seeds)

# segmentation
binary_segmentation = sitk.ConnectedThreshold(smoothed, seeds, lower=650, upper=1000)

# visualize result
vis.show_image_with_mask(input_image, binary_segmentation, 'Result')
