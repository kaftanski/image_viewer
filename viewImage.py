import argparse
import torch
import SimpleITK as sitk
import image_viewing
from GUI.viewer_utils import ImageMask

parser = argparse.ArgumentParser(description='3D Slice Image Viewer.')
parser.add_argument('img', type=str, help='path to input image')
parser.add_argument('-l', type=str, help='path to label image')
parser.add_argument('-p', type=str, help='path to landmarks (torch tensor object)')
args = parser.parse_args()

img = sitk.ReadImage(args.img)
if args.l:
    lbl = sitk.ReadImage(args.l)
    lbl = ImageMask(lbl)
else:
    lbl = None

if args.p:
    pts = torch.load(args.p, map_location=torch.device('cpu'))
    if pts.shape[0] == 3:
        pts = pts.transpose(0, 1)
    pts = (pts + 1) / 2 * torch.tensor(img.GetSize())
    pts = pts.tolist()
else:
    pts = None

image_viewing._show_image(img, window_title=args.img, blocking=True, mask=lbl, return_markers=False, points=pts)
